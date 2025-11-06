rm(list = ls())
.libPaths("/home/lblackma/R/x86_64-pc-linux-gnu-library/4.5/")
options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages("tidyverse", "stringr")
#remove.packages("INLAloggrowth")
remotes::install_github("ljb34/INLAlogisticgrowth", force = F)

library(inlabru)
library(Matrix)
library(sf)
library(dplyr)
library(INLAloggrowth)



out.lgcp <- simulate_loggrowth(growth = 0.8, carry.cap = 500, movement = 0.25, sigma = 0.1,
                               initial.pop = 250, initial.range = 0.3, initial.sigma=0.01, 
                               timesteps = 4,sample.type = "LGCP", boundaries = c(0,1), max.edge = 0.1, debug = T)

saveRDS(out.lgcp, "/home/lblackma/lgcp_data_gamma.RData")

library(dplyr)
library(tidyr)
print("Simulated with initial pop 250 and K = 500, observed:")
print(out.lgcp$animal_obs %>% group_by(time) %>% summarise(n = n()))

print("Mean of field on log and exp scale:")
print(out.lgcp$animal %>% group_by(time) %>% summarise(mean = mean(field, na.rm = T), em = mean(exp(field), na.rm = T)))


#fit initial year
bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
mesh_obs <- fm_mesh_2d(boundary = bnd,
                       max.edge = c(0.075,0.3), offset = c(-0.1,0.5))
bnd <- st_as_sf(bnd)
mesh_time <- fm_mesh_1d(loc = 1:4)
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(1, 0.1),
                              prior.range = c(0.01, 0.1))
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
initial.precision <- Q#Diagonal(mesh_obs$n,
#1/((fit0$summary.fixed$sd**2)+(fit0$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2)))+Q
#fit other years
priors <- list(cc = c(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,])/25,1/25),
               growth = c(log(0.75),0.15),move = c(0.25,0.1),sigma = c(log(0.1),0.1))
initial.linpoint <- log(logit.nest(exp(fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n]),
                                   0.8,nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,]),4)$x)


iterated.fit.lgcp <- iterate.fit.lgcp(data = out.lgcp$animal_obs, smesh = mesh_obs, tmesh = mesh_time,
                                      samplers = bnd,prior.mean = rep(fit0$summary.fixed$mean,mesh_obs$n)+fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                      prior.precision = initial.precision, priors = priors,
                                      max.iter = 10,gamma = 0.5,
                                      stop.crit = 0.01,initial.growth = log(0.75),
                                      initial.log.sigma = log(0.1),
                                      initial.move.const = 0.25,
                                      initial.linpoint = initial.linpoint, initial.carry.cap = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,])),
                                      options = list(verbose = T, control.inla = list(int.strategy = 'eb')), saveall = T)
iterated.fit.lgcp$data <- out.lgcp
iterated.fit.lgcp$initial <- fit0
saveRDS(iterated.fit.lgcp, "/home/lblackma/full_slurm_standard_gamma.RData")
pred.pixels <- fm_pixels(mesh_obs, mask = bnd, format = "sf")
pred.pixels.time <- fm_cprod(pred.pixels, data.frame(time = c(1:4)))
preds <- predict(iterated.fit.lgcp$fit, pred.pixels.time,
                 ~data.frame(loglambda = loggrow,
                             lambda = exp(loggrow)),
                 n.samples = 100)
saveRDS(preds, "/home/lblackma/full_slurm_standard_preds_gamma.RData")
