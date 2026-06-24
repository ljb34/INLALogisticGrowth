simulate_loggrowth_vary2 <- function(growth0, growth1, carry.cap0,carry.cap1, movement0, movement1, sigma,
                                    cov.range, cov.sigma,
                                    initial.pop,initial.range, initial.sigma, 
                                    timesteps, npoints = NULL, obs.sd=NULL,
                                    obs.prob = NULL, same.cov = T,
                                    sample.type = "LGCP", ncores = 1,
                                    boundaries = c(0,1), debug = F,
                                    max.edge = 0.05, nsurv = 3){
 
  #set up for simulation
  bnd_extended <- sf::st_as_sf(inlabru::spoly(data.frame(easting = c(boundaries[1], boundaries[2],boundaries[2],boundaries[1]), 
                                                         northing = c(boundaries[1], boundaries[1],boundaries[2],boundaries[2]))))
  hex_points <- fm_hexagon_lattice(bnd = bnd_extended, edge_len = 0.9*max.edge)
  smesh <- fmesher::fm_mesh_2d_inla(loc = hex_points, boundary = bnd_extended,
                                    max.edge = c(max.edge*1.1, 2*max.edge),
                                    offset = c(-0.01, (boundaries[2]-boundaries[1])))
  tmesh <- fmesher::fm_mesh_1d(loc = 0:timesteps)
  step.size <- 1
  if(debug) print("set up finished, generating first year")
  matern <-
    inla.spde2.pcmatern(smesh,
                        prior.sigma = c(0.1, 0.1),
                        prior.range = c(0.1, 0.1))
  initial_Q <- inla.spde.precision(matern,
                                   theta = log(c(initial.range, initial.sigma))) 
  prior.mean <- log(initial.pop) + inla.qsample(1, initial_Q)[,1]
  
  #components needed for model
  cov_Q <- inla.spde.precision(matern,theta = log(c(cov.range, cov.sigma)))
  if(same.cov){
    covariates <- data.frame(growth = inla.qsample(1, cov_Q)[, 1]) %>% 
      dplyr::mutate(carry.cap = growth, movement = growth)
  } else{
    covariates <- data.frame(growth = inla.qsample(1, cov_Q)[, 1],
                             carry.cap = inla.qsample(1, cov_Q)[, 1],
                             movement = inla.qsample(1, cov_Q)[, 1])
  }
  cov.grid <- sf::st_as_sf(expand.grid(
    easting = seq(boundaries[1],boundaries[2], by = 0.01),
    northing = seq(boundaries[1],boundaries[2], by = 0.01)), coords = c("easting", "northing"))
  cov.grid$growth <- fmesher::fm_evaluate(
    smesh,
    loc = cov.grid,
    field = covariates$growth)
  cov.grid$carry.cap <- fmesher::fm_evaluate(
    smesh,
    loc = cov.grid,
    field = covariates$carry.cap)
  cov.grid$movement <- fmesher::fm_evaluate(
    smesh,
    loc = cov.grid,
    field = covariates$movement)
  
  growth <- rep(exp(growth0 + growth1*covariates$growth),timesteps+1)
  carry.cap <- rep(exp(carry.cap0 + carry.cap1*covariates$carry.cap), timesteps+1)
  move.const <- rep(movement0 + movement1*covariates$movement, timesteps+1)
  print(summary(carry.cap))
  print(summary(growth))
  par = list(growth = growth, carry.cap = carry.cap, move.const = move.const, sigma = sigma)
  linpoint <- log(logit.nest(exp(prior.mean), growth[1:smesh$n], carry.cap[1:smesh$n], tmesh$n)$x)
  grad <- gradient_of_linpoint(linpoint, smesh, tmesh)#
  prior.precision <- initial_Q
  
  
  browser()
  cgen <- define.varying.cgeneric.loggrow.model(linpoint = linpoint, smesh = smesh,
                                                tmesh = tmesh, step.size = 1,
                                                prior.mean = prior.mean, 
                                                prior.precision = prior.precision,
                                                growth.formula = ~Intercept + growth_cov,
                                                carry.formula = ~Intercept + carry_cov,
                                                move.formula = ~Intercept + move_cov,
                                                growth_cov = covariates$growth,
                                                carry_cov = covariates$carry.cap,
                                                move_cov = covariates$movement,
                                                priors = NULL, grad = grad)
  
  theta = c(growth0, growth1, carry.cap0,carry.cap1, movement0, movement1, sigma)
  
  Q_mat <- INLAtools::cgeneric_Q(cgen, theta)
  
  
}

library(INLA)
library(inlabru)
library(INLAloggrowth)
library(fmesher)
library(sf)
library(dplyr)
library(ggplot2)

simulate_loggrowth_vary2(log(0.8), 0.2, log(1000), -0.5, 0.15, 0, 0.05,cov.range = 0.2,cov.sigma= 0.5, initial.pop = 500, 0.1, 0.1, 4,
                         max.edge = 0.5)

