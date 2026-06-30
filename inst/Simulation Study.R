library(INLA)
library(inlabru)
library(INLAloggrowth)
library(INLAspacetime)
library(sf)
library(dplyr)
#out.lgcp <- simulate_loggrowth(growth = 0.8, carry.cap = 1000, movement = 0.15, sigma = 5,
#                               initial.pop = 400, initial.range = 0.15, initial.sigma=0.05,
#                               timesteps = 10,sample.type = "LGCP", boundaries = c(0,1), debug = T,
#                               max.edge = 0.11)
out.lgcp <- readRDS("simstudy_data_gof5.RData")

dataobs <- filter(out.lgcp$animal_obs, time <=5)
#print(out.lgcp$animal_obs %>% group_by(time) %>% summarise(n = n()))

bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
bnd <- st_as_sf(bnd)
hex_points <- fm_hexagon_lattice(bnd = bnd, edge_len = 0.08)
mesh_obs <- fm_mesh_2d(locs = hex_points, boundary = bnd,
                       max.edge = c(0.05, 0.25), offset = c(-0.1, 1.5))

mesh_time <- fm_mesh_1d(loc = 1:5)
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(0.5, 0.05),
                              prior.range = c(0.1, 0.02))
cmp <- geometry ~ smooth(geometry, model = matern) +
  initial(1,model = "linear")-1

subdiv <- fm_subdivide(mesh_obs,1)
fit0 <- bru(cmp, out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,],domain = list(geometry = subdiv),
            family = "cp",samplers = bnd)

index <- min(which(stringr::str_sub(rownames(fit0$summary.fitted.values),8,8)!= "A"))
n.nodes <- fit0$misc$configs$nconfig
nodes <- data.frame(log.prob=rep(NA,n.nodes))
mat_list <- list()

for(i in 1:n.nodes){
  nodes[i,]<- fit0$misc$configs$config[[i]]$log.posterior
  Q <- fit0$misc$configs$config[[i]]$Q[1:mesh_obs$n, 1:mesh_obs$n]
  dQ <- diag(Q)
  Q <- Q + Matrix::t(Q)
  diag(Q) <- dQ
  mat_list[[i]] <- Q
}
nodes <- dplyr::mutate(nodes, weight = exp(log.prob)) %>%
  dplyr::mutate(weight.prob = weight/sum(weight))
Q <- Reduce("+", Map(function(m,w) w*m, mat_list,nodes$weight.prob))

initial.precision <-Diagonal(mesh_obs$n,
                             1/((fit0$summary.fixed$sd**2)+(fit0$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2)))%*%Q
priors <- list(cc = c(log(1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,])),0.25),
               growth = c(log(0.8),0.4),move = c(log(0.15),0.4),sigma = c(log(5),0.4))


summary(exp(fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n]))
mesh_locs <- st_as_sf(data.frame(x = mesh_obs$loc[,1], y = mesh_obs$loc[,2]), coords = c("x","y"))
pred0 <- predict(fit0, mesh_locs, ~initial+smooth)

initial.linpoint <- pred0$mean


for(i in 2:5){
 print(paste("Fitting year", i))
 cmp <- geometry ~ smooth(geometry, model = matern) +
   initial(1,model = "linear")-1
 fiti <- bru(cmp, out.lgcp$animal_obs[out.lgcp$animal_obs$time == i,],
             domain = list(geometry = subdiv),
             family = "cp",samplers = bnd)
 initial_fits[[i]] <- fiti
 predi <- predict(fiti, mesh_locs, ~initial+smooth)
 print(summary(predi))
 initial.linpoint <- c(initial.linpoint,predi$mean)
 rm(fiti)
 rm(cmp)
 rm(predi)
}


plot(exp(initial.linpoint))
print(summary(exp(initial.linpoint)))

initial.linpoint[initial.linpoint >= log(1500)] <- max(initial.linpoint[initial.linpoint < log(1500)], na.rm = T)
plot(exp(initial.linpoint))
print(summary(exp(initial.linpoint)))
print(summary(dataobs))


iterated.fit.lgcp <- iterate.cgeneric.fit.lgcp(data = dataobs, smesh = mesh_obs, tmesh = mesh_time,
                                              samplers = bnd,prior.mean = fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                              prior.precision = initial.precision, priors = priors,
                                              max.iter = 5,gamma = 0.25,
                                              stop.crit = 0.01, domain = list(geometry = subdiv, time = 1:5),
                                              #initial.linpoint = initial.linpoint,
                                              initial.growth = log(0.8),
                                              initial.log.sigma = log(5),
                                              initial.move.const = log(0.15),
                                              initial.carry.cap = log(1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,])),
                                              saveall = F, options = list(verbose = T, control.compute = list(dic = TRUE)), update.rule = 2)


#saveRDS(iterated.fit.lgcp, "/mnt/shared/scratch/lblackma/simstudy_loggrow_gof5_fine.RData")

pred.pixels <- fm_pixels(mesh_obs, mask = bnd, format = "sf")
pred.pixels.time <- fm_cprod(pred.pixels, data.frame(time = c(1:10)))
preds <- predict(iterated.fit.lgcp$fit, pred.pixels.time,
                ~data.frame(loglambda = loggrow,
                            lambda = exp(loggrow)),
                n.samples = 100)
#saveRDS(preds, "/mnt/shared/scratch/lblackma/simstudy_preds_loggrow_gof5_fine.RData")


######Logistic Growth functions ####
logit.growth <- function(x,r,k) {
  x[x<0] <- 0
  xnew <- x*k*exp(r)/(k-x + x*exp(r))
  return(xnew)
}

logit.nest.lik <- function(x0,r,k,n){
  x <- x0
  x <- x0*k*exp(r*n)/(k-x0 + x*exp(r*n))
  return(log(x))
}
######Diffusion ########
stmodel <- stModel.define(smesh = mesh_obs, ## spatial mesh
                          tmesh = mesh_time, ## temporal mesh
                          model = '121', ## model, see the paper
                          control.priors = list(
                            prs = c(0.05, 0.1), ## P(spatial range < 0.01) = 0.1
                            prt = c(1, 0.1), psigma = c(1, 0.1) ## P(sigma > 1) = 0.1
                          )
)
diffusion_cmp <- geometry + time ~ smooth(list(space = geometry, 
                                               time = time),
                                          model = stmodel) + 
  r(1, model = "linear",mean.linear = 0.8, prec.linear = 16)+
  k(1, model = "linear", mean.linear = 1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,]), prec.linear = (1/200)**2)+
  init0(1, model = "linear", mean.linear = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,])))


diffusion_lik1 <- bru_obs(formula = geometry+time ~init0 + smooth-1,
                          family = "cp",
                          data = dataobs[dataobs$time == 1,],
                          domain = list(geometry = mesh_obs, time = 1),
                          samplers = bnd)


diffusion_lik2 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                                  r, k, 1) + smooth-1,
                          family = "cp",
                          data = dataobs[dataobs$time == 2,],
                          domain = list(geometry = mesh_obs, time = 2),
                          samplers = bnd)
diffusion_lik3 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                                  r, k, 2) + smooth-1,
                          family = "cp",
                          data = dataobs[dataobs$time == 3,],
                          domain = list(geometry = mesh_obs, time = 3),
                          samplers = bnd)
diffusion_lik4 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                                  r, k, 3) + smooth-1,
                          family = "cp",
                          data = dataobs[dataobs$time == 4,],
                          domain = list(geometry = mesh_obs, time = 4),
                          samplers = bnd)
diffusion_lik5 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                                  r, k, 4) + smooth-1,
                          family = "cp",
                          data = dataobs[dataobs$time == 5,],
                          domain = list(geometry = mesh_obs, time = 5),
                          samplers = bnd)
diffusion_fit <- bru(diffusion_cmp, diffusion_lik1, diffusion_lik2,diffusion_lik3, diffusion_lik4,diffusion_lik5,
                     options = list(control.compute = list(dic = TRUE, waic = T),
                                    bru_initial = list(r = 0.8, k=1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,])),
                                    bru_max_iter = 50))
#saveRDS(diffusion_fit,file = "simstudy_diffusion_fit_gof5_prior1_fine.RData")
diffusion_preds <- predict(diffusion_fit, pred.pixels.time,
                           ~data.frame(
                             loglamb1 = init0 + smooth,
                             loglamb2 = logit.nest.lik(exp(init0),
                                                       r, k,1)+smooth,
                             loglamb3 = logit.nest.lik(exp(init0),
                                                       r, k, 2)+smooth,
                             loglamb4 = logit.nest.lik(exp(init0),
                                                       r,k, 3)+smooth,
                             loglamb5 = logit.nest.lik(exp(init0),
                                                       r, k, 4)+smooth,
                             loglamb6 = logit.nest.lik(exp(init0),
                                                       r, k, 5)+smooth,
                             loglamb7 = logit.nest.lik(exp(init0),
                                                       r, k, 6)+smooth,
                             loglamb8 = logit.nest.lik(exp(init0),
                                                       r, k, 7)+smooth,
                             loglamb9 = logit.nest.lik(exp(init0),
                                                       r, k, 8)+smooth,
                             loglamb10 = logit.nest.lik(exp(init0),
                                                        r, k, 9)+smooth),
                           
                           
                           n.samples = 100)
#saveRDS(diffusion_preds, "simstudy_diffusion_preds_gof5_prior1_fine.RData")
rm(diffusion_fit, diffusion_preds)
gc()

########## IID ###########
iid_cmp <- geometry + time ~ smooth(geometry, model = matern, group = time, 
                                    ngroup = 10,
                                    control.group = list(model = "iid")) + 
  r(1, model = "linear",mean.linear = 0.8, prec.linear = 16)+
  k(1, model = "linear", mean.linear = 1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,]), prec.linear = (1/200)**2)+
  init0(1, model = "linear", mean.linear = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,])))

iid_lik1 <- bru_obs(formula = geometry+time ~init0 + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 1,],
                    domain = list(geometry = mesh_obs, time = 1),
                    samplers = bnd)


iid_lik2 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                            r, k, 1) + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 2,],
                    domain = list(geometry = mesh_obs, time = 2),
                    samplers = bnd)
iid_lik3 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                            r, k, 2) + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 3,],
                    domain = list(geometry = mesh_obs, time = 3),
                    samplers = bnd)
iid_lik4 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                            r, k, 3) + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 4,],
                    domain = list(geometry = mesh_obs, time = 4),
                    samplers = bnd)

iid_lik5 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                            r, k, 4) + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 5,],
                    domain = list(geometry = mesh_obs, time = 5),
                    samplers = bnd)
iid_fit <- bru(iid_cmp, iid_lik1, iid_lik2,iid_lik3, iid_lik4, iid_lik5,
               options = list(control.compute = list(dic = TRUE, waic = T),
                              bru_initial = list(r = 0.8, k=1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,])),
                              bru_max_iter = 50))
#saveRDS(iid_fit,file = "simstudy_iid_fit_gof5_prior1_fine.RData")
iid_preds <- predict(iid_fit, pred.pixels.time,
                     ~data.frame(
                       loglamb1 = init0 + smooth,
                       loglamb2 = logit.nest.lik(exp(init0),
                                                 r, k,1)+smooth,
                       loglamb3 = logit.nest.lik(exp(init0),
                                                 r, k, 2)+smooth,
                       loglamb4 = logit.nest.lik(exp(init0),
                                                 r,k, 3)+smooth,
                       loglamb5 = logit.nest.lik(exp(init0),
                                                 r, k, 4)+smooth,
                       loglamb6 = logit.nest.lik(exp(init0),
                                                 r, k, 5)+smooth,
                       loglamb7 = logit.nest.lik(exp(init0),
                                                 r, k, 6)+smooth,
                       loglamb8 = logit.nest.lik(exp(init0),
                                                 r, k, 7)+smooth,
                       loglamb9 = logit.nest.lik(exp(init0),
                                                 r, k, 8)+smooth,
                       loglamb10 = logit.nest.lik(exp(init0),
                                                  r, k, 9)+smooth),
                     
                     
                     n.samples = 100)
#saveRDS(iid_preds, "simstudy_iid_preds_gof5_prior1_fine.RData")


############## AR1 ###############
ar1_cmp <- geometry + time ~ smooth(geometry, model = matern, group = time, 
                                    ngroup = 10,
                                    control.group = list(model = "ar1")) + 
  r(1, model = "linear",mean.linear = 0.8, prec.linear = 16)+
  k(1, model = "linear", mean.linear = 1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,]), prec.linear = (1/200)**2)+
  init0(1, model = "linear", mean.linear = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,])))

ar1_lik1 <- bru_obs(formula = geometry+time ~init0 + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 1,],
                    domain = list(geometry = mesh_obs, time = 1),
                    samplers = bnd)


ar1_lik2 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                            r, k, 1) + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 2,],
                    domain = list(geometry = mesh_obs, time = 2),
                    samplers = bnd)
ar1_lik3 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                            r, k, 2) + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 3,],
                    domain = list(geometry = mesh_obs, time = 3),
                    samplers = bnd)
ar1_lik4 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                            r, k, 3) + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 4,],
                    domain = list(geometry = mesh_obs, time = 4),
                    samplers = bnd)
ar1_lik5 <- bru_obs(formula = geometry+time ~logit.nest.lik(exp(init0),
                                                            r, k, 4) + smooth-1,
                    family = "cp",
                    data = dataobs[dataobs$time == 5,],
                    domain = list(geometry = mesh_obs, time = 5),
                    samplers = bnd)
ar1_fit <- bru(ar1_cmp, ar1_lik1, ar1_lik2,ar1_lik3,ar1_lik4,ar1_lik5,
               options = list(control.compute = list(dic = TRUE, waic = T),
                              bru_initial = list(r = 0.8, k=1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,])),
                              bru_max_iter = 50))
#saveRDS(ar1_fit,file = "simstudy_ar1_fit_gof5_prior1_fine.RData")
ar1_preds <- predict(ar1_fit, pred.pixels.time,
                     ~data.frame(
                       loglamb1 = init0 + smooth,
                       loglamb2 = logit.nest.lik(exp(init0),
                                                 r, k,1)+smooth,
                       loglamb3 = logit.nest.lik(exp(init0),
                                                 r, k, 2)+smooth,
                       loglamb4 = logit.nest.lik(exp(init0),
                                                 r,k, 3)+smooth,
                       loglamb5 = logit.nest.lik(exp(init0),
                                                 r, k, 4)+smooth,
                       loglamb6 = logit.nest.lik(exp(init0),
                                                 r, k, 5)+smooth,
                       loglamb7 = logit.nest.lik(exp(init0),
                                                 r, k, 6)+smooth,
                       loglamb8 = logit.nest.lik(exp(init0),
                                                 r, k, 7)+smooth,
                       loglamb9 = logit.nest.lik(exp(init0),
                                                 r, k, 8)+smooth,
                       loglamb10 = logit.nest.lik(exp(init0),
                                                  r, k, 9)+smooth),
                     
                     
                     n.samples = 100)
#saveRDS(ar1_preds, "simstudy_ar1_preds_gof5_prior1_fine.RData")
