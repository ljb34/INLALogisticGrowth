#.libPaths("mnt/apps/users/lblackma/RPackages/")
#install.packages(c("remotes", "sf", "dplyr", "tidyr", "stringr"))
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
rm(list = ls())
options(repos = c(CRAN = "https://cloud.r-project.org"))
library(INLA)
remove.packages("INLAloggrowth")
remotes::install_github("ljb34/INLAlogisticgrowth@cgeneric-blocks", force = T)
library(inlabru)
library(INLAloggrowth)
library(sf)
library(dplyr)
out.lgcp <- simulate_loggrowth(growth = 0.8, carry.cap = 500, movement = 0.25, sigma = 0.2,
                               initial.pop = 250, initial.range = 0.3, initial.sigma=0.01, 
                               timesteps = 4,sample.type = "LGCP", boundaries = c(0,1), debug = T,
                               max.edge = 0.1)
saveRDS(out.lgcp, "cgensmalldata.RData")
#out.lgcp <- readRDS("/mnt/shared/scratch/lblackma/cgensmalldata.RData")
#fit initial year
if(nrow(out.lgcp$animal_obs) <= 200){
  print("Something's gone wrong, retrying")
  out.lgcp <- simulate_loggrowth(growth = 0.8, carry.cap = 500, movement = 0.25, sigma = 0.1,
                               initial.pop = 250, initial.range = 0.3, initial.sigma=0.025, 
                               timesteps = 4,sample.type = "LGCP", boundaries = c(0,1), debug = T,
                               max.edge = 0.1)
  saveRDS(out.lgcp, "/mnt/shared/scratch/lblackma/cgensmalldata.RData")
  if(nrow(out.lgcp$animal_obs) <= 200){
    print("Failed to simulate data")
    quit("no")
    }
}
print(nrow(out.lgcp$animal_obs))
bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
mesh_obs <- fm_mesh_2d(boundary = bnd,
                       max.edge = c(0.075,0.5), offset = c(-0.1, 0.75))
bnd <- st_as_sf(bnd)
mesh_time <- fm_mesh_1d(loc = 1:4)
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(1, 0.1),
                              prior.range = c(0.025, 0.1))
cmp <- geometry ~ smooth(geometry, model = matern) +
  initial(1,model = "linear", mean.linear = log(length(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,])))-1

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

initial.precision <-Diagonal(mesh_obs$n,
                             1/((fit0$summary.fixed$sd**2)+(fit0$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2)))+Q
#fit other years
priors <- list(cc = c(log(450),0.2),
               growth = c(log(0.8),0.2),move = c(0.25,0.1),sigma = c(log(0.1),0.25))
print(priors)
summary(exp(fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n]))
initial.linpoint <- log(logit.nest(exp(fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n]), 0.8, 500,4)$x)
summary(exp(initial.linpoint))
iterated.fit.lgcp <- iterate.cgeneric.fit.lgcp(data = out.lgcp$animal_obs, smesh = mesh_obs, tmesh = mesh_time,
                                               samplers = bnd,prior.mean = fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                               prior.precision = initial.precision, priors = priors,
                                               max.iter = 20,gamma = 0.5,
                                               stop.crit = 0.1,
                                               initial.linpoint = NULL, initial.growth = log(0.8),
                                               initial.log.sigma = log(0.1),
                                               initial.move.const = 0.25,
                                               initial.carry.cap = log(500),
                                               saveall = T, options = list(verbose = T, control.inla = list(int.strategy = 'eb')))
iterated.fit.lgcp$data <- out.lgcp
iterated.fit.lgcp$initial <- fit0
saveRDS(iterated.fit.lgcp, "/mnt/shared/scratch/lblackma/cgen_small_blocks.RData")
pred.pixels <- fm_pixels(mesh_obs, mask = bnd, format = "sf")
pred.pixels.time <- fm_cprod(pred.pixels, data.frame(time = c(1:4)))
preds <- predict(iterated.fit.lgcp$fit, pred.pixels.time,
                 ~data.frame(loglambda = loggrow,
                             lambda = exp(loggrow)),
                 n.samples = 100)
saveRDS(preds, "/mnt/shared/scratch/lblackma/cgen_small_preds_blocks.RData")


