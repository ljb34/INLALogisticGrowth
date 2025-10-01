#' Cgeneric unit test
library(INLAloggrowth)
library(sf)
out.lgcp <- simulate_loggrowth(growth = 1, k = 500, movement = 1, sigma = 20,
                               initial = 250,timesteps = 4,sample.type = "LGCP")

#fit initial year
bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
mesh_obs <- fm_mesh_2d(boundary = bnd,
                       max.edge = c(0.05,1))
bnd <- st_as_sf(bnd)
mesh_time <- fm_mesh_1d(loc = 1:4)
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(0.1, 0.1),
                              prior.range = c(0.1, 0.1))
cmp <- geometry ~ smooth(geometry, model = matern) +
  initial(1,model = "linear", mean.linear = log(length(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,])))-1

fit0 <- bru(cmp, out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,],domain = list(geometry = mesh_obs),
            family = "cp",samplers = bnd)
#Find fitted values on mesh points
library(stringr)
index <- min(which(str_sub(rownames(fit0$summary.fitted.values),8,8)!= "A"))

initial.variance <- Diagonal(mesh_obs$n, (fit0$summary.fixed$sd**2)+(fit0$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2))
#fit other years
priors <- list(cc = c(log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,])),1),
               growth = c(1,1),move = c(1,1),sigma = c(log(20),1))

iterated.fit.lgcp <- iterate.cgeneric.fit.lgcp(data = out.lgcp$animal_obs, smesh = mesh_obs, tmesh = mesh_time,
                                               samplers = bnd,prior.mean = fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                               prior.variance = initial.variance, priors = priors,
                                               max.iter = 100,gamma = 0.5,
                                               stop.crit = 0.01,
                                               initial.linpoint = NULL, initial.growth = 0.8,
                                               initial.carry.cap = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,])),
                                               verbose = T)
iterated.fit.lgcp$data <- out.lgcp
saveRDS(iterated.fit.lgcp, "LogGrowth/cgeneric.RData")
pred.pixels <- fm_pixels(mesh_obs, mask = bnd, format = "sf")
pred.pixels.time <- fm_cprod(pred.pixels, data.frame(time = c(1:4)))
preds <- predict(iterated.fit.lgcp$fit, pred.pixels.time,
                 ~data.frame(time = 1:4, loglambda = loggrow,
                             lambda = exp(loggrow)),
                 n.samples = 100)
saveRDS(preds, "LogGrowth/cgeneric_preds.RData")


