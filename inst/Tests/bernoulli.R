rm(list = ls())
.libPaths("/home/lblackma/R/x86_64-pc-linux-gnu-library/4.5/")
options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages("tidyverse", "stringr")
remove.packages("INLAloggrowth")
remotes::install_github("ljb34/INLAlogisticgrowth", force = F)

library(inlabru)
library(Matrix)
library(sf)
library(dplyr)
library(INLAloggrowth)
bernoulli.obs <- simulate_loggrowth(growth = 0.8, carry.cap = 10, movement = 0.25, sigma = 0.1,
                                    initial.pop = 1, initial.range = 0.3, initial.sigma=0.01, 
                                    timesteps = 4,sample.type = "Bernoulli", npoints = 100, obs.prob = 0.075,
                                    boundaries = c(0,1),max.edge = 0.1, nsurv = 1, debug = T)

#fit initial year
bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
mesh_obs <- fm_mesh_2d(boundary = bnd,
                       max.edge = c(0.05,0.4), offset = c(-0.1,0.6))
bnd <- st_as_sf(bnd)
mesh_time <- fm_mesh_1d(loc = 1:4)
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(1, 0.1),
                              prior.range = c(0.025, 0.1))
cmp <- obs ~ smooth(geometry, model = matern) +
  initial(1,model = "linear")-1

fit0 <- bru(cmp, bernoulli.obs$animal_obs[bernoulli.obs$animal_obs$time == 1,],domain = list(geometry = mesh_obs),
            family = "binomial", Ntrials = 1)
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
priors <- list(cc = c(90,0.15),
               growth = c(log(0.7),0.15),move = c(0.25,0.1),sigma = c(log(0.1),0.1))
initial.linpoint <- log(logit.nest(exp(plogis(fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n])),
                                   0.8,10,4)$x)

iterate.bernoulli.fit <- iterate.fit.custom(obs ~ -1 + prob(1, model = "linear", mean.linear = 0.05), 
                                            data = bernoulli.obs$animal_obs,
                                            family = "binomial", smesh = mesh_obs,
                                            tmesh = mesh_time, samplers = bnd,
                                            prior.mean = exp(plogis(fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n])),
                                            prior.precision = initial.precision,
                                            max.iter = 10, gamma = 0.25, priors = priors,
                                            initial.growth = log(0.7), initial.carry.cap = log(10),
                                            initial.move.const = 0.25,
                                            initial.log.sigma = log(0.1), options = list(verbose = T), debug = T,
                                            method = "rgeneric")

iterate.bernoulli.fit$data <- bernoulli.obs
iterate.bernoulli.fit$initial <- fit0
saveRDS(iterate.bernoulli.fit, "/home/lblackma/bernoulli_fit.RData")
pred.pixels <- fm_pixels(mesh_obs, mask = bnd, format = "sf")
pred.pixels.time <- fm_cprod(pred.pixels, data.frame(time = c(1:4)))
preds <- predict(iterate.bernoulli.fit$fit, pred.pixels.time,
                 ~data.frame(loglambda = loggrow,
                             lambda = exp(loggrow)),
                 n.samples = 100)
saveRDS(preds, "/home/lblackma/bernoulli_preds.RData")