remotes::install_github("ljb34/INLAlogrowth", ref = "cgeneric")
library(INLAloggrowth)
library(dplyr)
library(INLA)
library(fmesher)
library(sf)
out.lgcp <- simulate_loggrowth(growth = 0.8, k = 500, movement = 0.5, sigma = 20,
                               initial = 250,timesteps = 4,sample.type = "LGCP")

#fit initial year
bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
mesh_obs <- fm_mesh_2d_inla(boundary = bnd,
                            max.edge = c(0.07,1), offset = c(0,1))
bnd <- st_as_sf(bnd)
mesh_time <- fm_mesh_1d(loc = 1:4)
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(0.1, 0.1),
                              prior.range = c(0.1, 0.1))
cmp <- geometry ~ smooth(geometry, model = matern) +
  initial(1,model = "linear", mean.linear = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,])))-1

fit0 <- bru(cmp, out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,],domain = list(geometry = mesh_obs),
            family = "cp",samplers = bnd)
#Find fitted values on mesh points
library(stringr)
index <- min(which(str_sub(rownames(fit0$summary.fitted.values),8,8)!= "A"))
n.nodes <- fit0$misc$configs$nconfig
nodes <- data.frame(log.prob=rep(NA,n.nodes))
mat_list <- list()

for(i in 1:n.nodes){
  nodes[i,]<- fit0$misc$configs$config[[i]]$log.posterior
  mat_list[[i]] <- fit0$misc$configs$config[[i]]$Q[1:mesh_obs$n,1:mesh_obs$n]
}
nodes <- dplyr::mutate(nodes, weight = exp(log.prob)) %>%
  dplyr::mutate(weight.prob = weight/sum(weight))
Q <- Reduce("+", Map(function(m,w) w*m, mat_list,nodes$weight.prob))

initial.precision <- Diagonal(mesh_obs$n,
                              1/((fit0$summary.fixed$sd**2)+(fit0$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2)))+Q
#fit other years
priors <- list(cc = c(log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,])),1),
  growth = c(0.8,0.5),move = c(0.8,0.5),sigma = c(log(1),1))
initial.linpoint <- log(logit.nest(exp(fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n]),
                                   0.8,50,5)$x)
iterated.fit.lgcp <- iterate.fit.lgcp(data = out.lgcp$animal_obs, smesh = mesh_obs, tmesh = mesh_time,
                                      samplers = bnd,prior.mean = rep(fit0$summary.fixed$mean,mesh_obs$n)+fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                      prior.precision = initial.precision, priors = priors,
                                      max.iter = 100,gamma = 0.5,
                                      stop.crit = 0.01,initial.growth = 0.8,
                                      initial.carry.cap = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,])),
                                      initial.linpoint = initial.linpoint,
                                      verbose = T)
iterated.fit.lgcp$data <- out.lgcp
saveRDS(iterated.fit.lgcp, "LogGrowth/newsim.RData")
pred.pixels <- fm_pixels(mesh_obs, mask = bnd, format = "sf")
pred.pixels.time <- fm_cprod(pred.pixels, data.frame(time = c(1:5)))
preds <- predict(iterated.fit.lgcp$fit, pred.pixels.time,
                 ~data.frame(loglambda = loggrow,
                             lambda = exp(loggrow)),
                 n.samples = 100)
saveRDS(preds, "LogGrowth/newsim_preds.RData")