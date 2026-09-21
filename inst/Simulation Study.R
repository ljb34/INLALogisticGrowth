library(INLA)
library(inlabru)
library(INLAloggrowth)
library(INLAspacetime)
library(sf)
library(dplyr)
#Parameters used for simulating data. Not run as it takes ~5 minutes
#out.lgcp <- simulate_loggrowth(growth = 0.8, carry.cap = 1000, movement = 0.15, sigma = 5,
#                               initial.pop = 400, initial.range = 0.15, initial.sigma=0.05,
#                               timesteps = 10,sample.type = "LGCP", boundaries = c(0,1), debug = T,
#                               max.edge = 0.11)
out.lgcp <- readRDS("simstudy_data_gof5.RData")

dataobs <- filter(out.lgcp$animal_obs, time <=5)

#Set up for model fitting
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

#Fit first year to obtain prior
fit0 <- bru(cmp, out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,],domain = list(geometry = subdiv),
            family = "cp",samplers = bnd)

#Extract prior precision matrix
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

mesh_locs <- st_as_sf(data.frame(x = mesh_obs$loc[,1], y = mesh_obs$loc[,2]), coords = c("x","y"))

#Predict initial mean
pred0 <- predict(fit0, mesh_locs, ~initial+smooth)

initial.linpoint <- pred0$mean

#Fit years 2-5 to get starting guess for linearisation point
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
#Remove any extreme values for a smoother starting guess
initial.linpoint[initial.linpoint >= log(1500)] <- max(initial.linpoint[initial.linpoint < log(1500)], na.rm = T)

#Fit logistic growth model
iterated.fit.lgcp <- iterate.fit( formula = ~-1, #No intercept
                                  family = "cp",
                                  data = dataobs, smesh = mesh_obs, tmesh = mesh_time,
                                              samplers = bnd,prior.mean = fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                              prior.precision = initial.precision, priors = priors,
                                              max.iter = 5,gamma = 0.25,
                                              stop.crit = 0.01, domain = list(geometry = subdiv, time = 1:5),
                                              initial.linpoint = initial.linpoint,
                                              initial.growth = log(0.8),
                                              initial.log.sigma = log(5),
                                              initial.move.const = log(0.15),
                                              initial.carry.cap = log(1.25*nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 5,])),
                                              saveall = F, options = list(verbose = T, control.compute = list(dic = TRUE)), update.rule = 2)

#Generate Predictions
pred.pixels <- fm_pixels(mesh_obs, mask = bnd, format = "sf")
pred.pixels.time <- fm_cprod(pred.pixels, data.frame(time = c(1:10)))
preds <- predict(iterated.fit.lgcp$fit, pred.pixels.time,
                ~data.frame(loglambda = loggrow,
                            lambda = exp(loggrow)),
                n.samples = 100)



######Comparison with other methods -----------
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
stmodel <- stModel.define(smesh = mesh_obs, 
                          tmesh = mesh_time, 
                          model = '121', 
                          control.priors = list(
                            prs = c(0.05, 0.1), 
                            prt = c(1, 0.1), psigma = c(1, 0.1)))
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
##########Plots ##################
library(patchwork)
#Nice colour scale for plots
colsc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    limits = range(..., na.rm = TRUE)
  )
}
#Simulated Data
ggplot()+
  gg(out.lgcp$animal, aes(fill = exp(field)), geom = "tile")+
  gg(data$animal_obs)+
  facet_wrap(~time)+
  colsc(exp(out.lgcp$animal$field))+
  labs(fill = "Intensity")

#True intensity and predicted intensity for logistic growth model
ggplot()+
  gg(out.lgcp$animal[out.lgcp$animal$time <= 5,], aes(fill = exp(field)), geom = "tile")+
  facet_wrap(~time)+
  colsc(exp(out.lgcp$animal$field[out.lgcp$animal$time <= 5]), preds$lambda$median[preds$lambda$time <= 5])+
  ggtitle("True Intensity")+
  ggplot()+
  gg(preds$lambda[preds$lambda$time <= 5,], aes(fill = median), geom = "tile")+
  facet_wrap(~time)+
  colsc(exp(out.lgcp$animal$field[out.lgcp$animal$time <= 5]), preds$lambda$median[preds$lambda$time <= 5])+
  ggtitle("Predicted Intensity")

#Data processing to plot all models together
truth <- data$animal
truth$type <- "Truth"
preds$loglambda$type <- "Logistic"
for(i in 1:10){
  ar1preds[[i]]<- filter(ar1preds[[i]], time == i)
  ar1preds[[i]]$type <- "AR1"
  iidpreds[[i]]<- filter(iidpreds[[i]], time == i)
  iidpreds[[i]]$type <- "IID"
  diffpreds[[i]]<- filter(diffpreds[[i]], time == i)
  diffpreds[[i]]$type <- "Diffusion"
}
combpreds <- rbind(do.call(rbind, ar1preds), do.call(rbind, iidpreds), do.call(rbind, diffpreds),
                   preds$loglambda)
#Plot of truth and predictions for all models
ggplot()+
  gg(combpreds[combpreds$time <= 10,], aes(fill = exp(mean)), geom = "tile")+
  gg(truth[truth$time <=10,], aes(fill = exp(field)), geom = "tile")+
  facet_grid(rows = vars(factor(type, levels = c("Truth", "IID", "AR1", "Diffusion", "Logistic"))), 
             cols = vars(time),
             switch = "y")+
  colsc(exp(truth$field), exp(combpreds$mean[combpreds$time <=10]))+
  labs(fill = "Mean Intensity")+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

preds$lambda$type <- "Logistic"

#Plot intensity and uncertainty in same figure
intens <- ggplot()+
  gg(preds$lambda[preds$loglambda$time <= 5,], aes(fill = median), geom = "tile")+
  gg(truth[truth$time <=5,], aes(fill = exp(field)), geom = "tile")+
  facet_grid(rows = vars(factor(type, levels = c("Truth","Logistic"))), cols = vars(time), switch = "y")+
  colsc(exp(truth$field[truth$time <= 5]))+
  labs(fill = "Intensity")+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        strip.text.y = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
uncertainty <- preds$lambda
uncertainty$type <- "Uncertainty"
#uncertainty <- mutate(uncertainty, realsd = sqrt((exp(sd**2)-1)*exp(2*log(mean) + sd**2)))
sd <- ggplot()+
  gg(uncertainty[uncertainty$time <= 5,], aes(fill = sd), geom = "tile")+
  facet_grid(rows = vars(type), cols = vars(time), switch = "y")+
  scale_fill_continuous(type = "viridis")+
  #colsc(uncertainty$sd[uncertainty$time <= 5])+
  labs(fill = "SD")+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        strip.text.y = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

intens /sd

#Out of sample predictions
#Put on same scale
samescale <- c(max(preds$lambda$q0.975[preds$lambda$time > 5]), min(preds$lambda$q0.025[preds$lambda$time > 5]))
truth10 <- ggplot()+
  gg(truth[truth$time >5,], aes(fill = exp(field)), geom = "tile")+
  facet_grid(rows = vars(type), cols = vars(time), 
             switch = "y")+
  colsc(samescale)+
  labs(fill = "Intensity")+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

meanpreds <- preds$lambda %>% mutate(type = "Mean")
meanplot <-  ggplot()+
  gg(meanpreds[meanpreds$time > 5,], aes(fill = mean), geom = "tile")+
  facet_grid(rows = vars(type), cols = vars(time), 
             switch = "y")+
  colsc(samescale)+
  labs(fill = "Intensity")+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

q0.025preds <- preds$lambda %>% mutate(type = "q0.025")
q0.025plot <-  ggplot()+
  gg(q0.025preds[q0.025preds$time > 5,], aes(fill = q0.025), geom = "tile")+
  facet_grid(rows = vars(type), cols = vars(time), 
             switch = "y")+
  colsc(samescale)+
  labs(fill = "Intensity")+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

q0.975preds <- preds$lambda %>% mutate(type = "q0.975")
q0.975plot <-  ggplot()+
  gg(q0.975preds[q0.975preds$time > 5,], aes(fill = q0.975), geom = "tile")+
  facet_grid(rows = vars(type), cols = vars(time), 
             switch = "y")+
  colsc(samescale)+
  labs(fill = "Intensity")+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

truth10 / meanplot / q0.025plot / q0.975plot


###### Similarity scores #####
r <- 0.005

truth_pred <- combpreds %>%
  group_by(type, time) %>%
  group_modify(~{
    
    preds_t <- .x
    
    truth_t <- truth %>% select(-type) %>%
      filter(time == .y$time) %>%
      mutate(lambda = exp(field))
    
    nb <- st_is_within_distance(truth_t, preds_t,
                                dist = sqrt(2) * r + 1e-8)
    
    truth_t$pred_mean <- sapply(nb, function(i) {
      if (length(i) == 0) NA_real_
      else mean(exp(preds_t$mean[i]))
    })
    
    
    select(truth_t, -time)
  }) %>%
  ungroup() %>% st_as_sf()
dx <- 0.01
#Calculate SSIM
SSIM <- truth_pred %>%
  group_by(time, type) %>%
  group_modify(~{
    nb <- st_is_within_distance(.x, .x, dist = sqrt(2) * dx + 1e-8)
    
    # Remove self
    nb <- Map(setdiff, nb, seq_len(nrow(.x)))
    
    .x$neighbour_mean_pred <- sapply(nb, function(i) {
      if (length(i) == 0) NA_real_
      else mean(.x$pred_mean[i])
    })
    
    .x$sigma2_pred <- sapply(nb, function(i){
      if (length(i) == 0) NA_real_
      else var(.x$pred_mean[i])
    })
    
    .x$neighbour_mean_truth <- sapply(nb, function(i) {
      if (length(i) == 0) NA_real_
      else mean(.x$lambda[i])
    })
    
    .x$sigma2_truth <- sapply(nb, function(i){
      if (length(i) == 0) NA_real_
      else var(.x$lambda[i])
    })
    .x
  }) %>%
  ungroup() %>% st_as_sf()

SSIM <- SSIM %>%
  group_by(time, type) %>%
  group_modify(~{
    nb <- st_is_within_distance(.x, .x, dist = sqrt(2) * dx + 1e-8)
    
    # Remove self
    nb <- Map(setdiff, nb, seq_len(nrow(.x)))
    .x$sigma_ab <- sapply(nb, function(i){
      if (length(i) == 0) NA_real_
      else mean((.x$pred_mean-.x$neighbour_mean_pred)*(.x$lambda - .x$neighbour_mean_truth))
    })
    .x
  }) %>%
  ungroup() %>% st_as_sf()
R <- max(SSIM$lambda, SSIM$pred_mean)-min(SSIM$lambda, SSIM$pred_mean)

SSIM <- SSIM %>% group_by(time, type) %>%
  group_modify(~{
    .x$R <- max(.x$lambda, .x$pred_mean)-min(.x$lambda, .x$pred_mean)
    .x
  })
c1 = (0.01*R)**2
c2 = (0.03*R)**2
c3 = c2/2
SSIM <- SSIM %>% mutate(SIM = (2*neighbour_mean_pred*neighbour_mean_truth + (0.01*R)**2)/(neighbour_mean_pred**2 + neighbour_mean_truth**2 + (0.01*R)**2),
                        SIV = (2*sqrt(sigma2_pred*sigma2_truth)+(0.03*R)**2)/(sigma2_pred + sigma2_truth + (0.03*R)**2),
                        SIP = (sigma_ab + 0.5*(0.03*R)**2)/(sqrt(sigma2_pred*sigma2_truth)+0.5*(0.03*R)**2)) %>%
  mutate(SSIM = SIM*SIV*SIP)


SSIM %>% group_by(type) %>%summarise(meanSSIM = mean(SSIM), 
                                     meanSIM = mean(SIM), 
                                     meanSIV = mean(SIV),
                                     meanSIP = mean(SIP))
# Scores ------------------------------------------------------------------

# CRPS ---------------
#AR1
ar1_scorepreds <- predict(ar1fit, truth,
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
                          
                          
                          n.samples = 1000)


for(i in 1:10){
  ar1_scorepreds[[i]]<- filter(ar1_scorepreds[[i]], time == i)
}

ar1_score <- do.call(rbind, ar1_scorepreds)
maxK <- max(truth$field)+4*sqrt(max(truth$field))
minK <- min(truth$field)-4*sqrt(min(truth$field))

#check minKK and maxKK are ok
summary(pnorm(maxK,ar1_score$mean, ar1_score$sd)) #want to be close to 1
summary(pnorm(minK,ar1_score$mean, ar1_score$sd)) #want to be close to 0

kk <- seq(minK, maxK, length.out = 1000)



ar1_CRPS <- (pnorm(kk[1],ar1_score$mean, ar1_score$sd) -
               as.numeric(truth$field <= kk[1]))**2

for(k in kk){
  ar1_CRPS <- ar1_CRPS + (pnorm(k,ar1_score$mean, ar1_score$sd) -
                            as.numeric(truth$field <= k))**2
}
summary(ar1_CRPS)

###IID
iid_scorepreds <- predict(iidfit, truth,
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
                          
                          
                          n.samples = 1000)


for(i in 1:10){
  iid_scorepreds[[i]]<- filter(iid_scorepreds[[i]], time == i)
}

iid_score <- do.call(rbind, iid_scorepreds)

iid_CRPS <- (pnorm(kk[1],iid_score$mean, iid_score$sd) -
               as.numeric(truth$field <= kk[1]))**2

for(k in kk){
  iid_CRPS <- iid_CRPS + (pnorm(k,iid_score$mean, iid_score$sd) -
                            as.numeric(truth$field <= k))**2
}
summary(iid_CRPS)

#####Diffusion

diff_scorepreds <- predict(difffit, truth,
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
                           
                           
                           n.samples = 1000)


for(i in 1:10){
  diff_scorepreds[[i]]<- filter(diff_scorepreds[[i]], time == i)
}

diff_score <- do.call(rbind, diff_scorepreds)

diff_CRPS <- (pnorm(kk[1],diff_score$mean, diff_score$sd) -
                as.numeric(truth$field <= kk[1]))**2

for(k in kk){
  diff_CRPS <- diff_CRPS + (pnorm(k,diff_score$mean, diff_score$sd) -
                              as.numeric(truth$field <= k))**2
}
summary(diff_CRPS)

#Logistic growth
loggrow_score <- predict(iterated.fit.lgcp$fit, truth, ~loggrow)

loggrow_score$type <- "Logistic Growth"
loggrow_CRPS <- (pnorm(kk[1],loggrow_score$mean, loggrow_score$sd) -
                as.numeric(truth$field <= kk[1]))**2

for(k in kk){
  loggrow_CRPS <- loggrow_CRPS + (pnorm(k,loggrow_score$mean, loggrow_score$sd) -
                              as.numeric(truth$field <= k))**2
}

ar1_score$CRPS <- ar1_CRPS
iid_score$CRPS <- iid_CRPS
diff_score$CRPS <- diff_CRPS
loggrow_score$CRPS <- loggrow_CRPS
#MSE ------------
ar1_MSE <- (ar1_score$mean - truth$field)**2
summary(ar1_MSE)
iid_MSE <- (iid_score$mean - truth$field)**2
summary(iid_MSE)
diff_MSE <- (diff_score$mean - truth$field)**2
summary(diff_MSE)
loggrow_MSE <- (loggrow_score$mean - truth$field)**2
summary(loggrow_MSE)

ar1_score$MSE <- ar1_MSE
iid_score$MSE <- iid_MSE
diff_score$MSE <- diff_MSE
loggrow_score$MSE <- loggrow_MSE

ar1_MAE <- abs(ar1_score$median - ar1_score$field)
summary(ar1_MAE)
iid_MAE <- abs(iid_score$median - truth$field)
summary(iid_MAE)
diff_MAE <- abs(diff_score$median - truth$field)
summary(diff_MAE)
loggrow_MAE <- abs(loggrow_score$median - truth$field)
summary(loggrow_MAE)
ar1_score$MAE <- ar1_MAE
iid_score$MAE <- iid_MAE
diff_score$MAE <- diff_MAE
loggrow_score$MAE <- loggrow_MAE

summary(ar1_score)
summary(iid_score)
summary(diff_score)



ar1_score$type <- "AR1"
iid_score$type <- "IID"
diff_score$type <- "Diffusion"
scores <- rbind(ar1_score, iid_score, diff_score, loggrow_score)
scores %>% mutate(insample = time <= 5) %>% group_by(type, insample) %>% 
  summarise(meanCRPS = mean(CRPS), meanMSE = mean(MSE), meanMAE = mean(MAE))
