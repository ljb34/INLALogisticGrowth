#' Cgeneric unit test
library(INLAloggrowth)
library(sf)
out.lgcp <- simulate_loggrowth(growth = 1, k = 100, movement = 1, sigma = 1,
                               initial = 50,timesteps = 3,sample.type = "LGCP")

#fit initial year
bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
mesh_obs <- fm_mesh_2d(boundary = bnd,
                            max.edge = c(0.5,1), max.n = 20)
bnd <- st_as_sf(bnd)
mesh_time <- fm_mesh_1d(loc = 1:3)
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(0.1, 0.1),
                              prior.range = c(0.1, 0.1))
cmp <- geometry ~ smooth(geometry, model = matern) +
  initial(1,model = "linear", mean.linear = length(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,]))-1

fit0 <- bru(cmp, out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,],domain = list(geometry = mesh_obs),
            family = "cp",samplers = bnd)
#Find fitted values on mesh points
library(stringr)
index <- min(which(str_sub(rownames(fit0$summary.fitted.values),8,8)!= "A"))

initial.variance <- Diagonal(mesh_obs$n, (fit0$summary.fixed$sd**2)+(fit0$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2))
#fit other years
priors <- list(cc = c(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 3,]),50),
               growth = c(1,1),move = c(1,1),sigma = c(log(20),1))

iterated.fit.lgcp <- iterate.cgeneric.fit.lgcp(data = out.lgcp$animal_obs, smesh = mesh_obs, tmesh = mesh_time,
                                      samplers = bnd,prior.mean = fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                      prior.variance = initial.variance, priors = priors,
                                      max.iter = 100,gamma = 0.5,
                                      stop.crit = 0.01,
                                      initial.linpoint = NULL, initial.growth = 0.8,
                                      initial.carry.cap = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 3,])),
                                      verbose = T)


