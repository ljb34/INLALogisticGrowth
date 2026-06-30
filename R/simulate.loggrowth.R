#'Simulate Spatial Logistic Growth Data
#'
#'@description Simulates observations and underlying field for spatial logistic growth in a square area. 
#'
#'@param growth,k,movement,sigma parameters for spatial logistic growth model
#'@param initial starting population size
#'@param timesteps Number of years to simulate
#'@param sample.type One of "LGCP" or "Gaussian". Default LGCP
#'@param boundaries Lower and upper limit of each side of the square. Defaults to sides of length 1
#'@param npoints For \code{sample.type="Gaussian"} only. Number of points to sample.
#'@param obs.sd For \code{sample.type="Gaussian"} only. Standard deviation of observation process
#'@param ncores Optional. Number of cores to use if running in parallel. 
#'
#'
#'@returns List containing underlying field (animal), observations (animal_obs), and 
#'field defined on mesh nodes for debugging purposes (animal_field). 
#'
#'@examples test.data <- simulate.loggrowth(growth = 1,k = 200,movement = 1,
#'sigma = 10,initial = 75,timesteps = 3)
#'if(require("ggplot2")){
#'ggplot()+
#'gg(test.data$animal_obs)+facet_wrap(~time)+ggtitle("Observed animals")
#'ggplot()+
#'gg(test.data$animal, aes(fill = field), geom = "tile")+facet_wrap(~time)+
#'ggtitle("Underlying field")
#'}
#'
#'@export
simulate_loggrowth <- function(growth, carry.cap, movement, sigma, 
                               initial.pop,initial.range, initial.sigma, 
                               timesteps, npoints = NULL, obs.sd=NULL,
                               obs.prob = NULL,
                               sample.type = "LGCP", ncores = 1,
                               boundaries = c(0,1), debug = F,
                               max.edge = 0.05, nsurv = 3){
  #browser()
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
  
  if(debug){
    print("Defining model")
  }
  #components needed for model
  theta <- c(log(growth), log(carry.cap), log(movement), log(sigma))
  linpoint <- log(logit.nest(exp(prior.mean), growth, carry.cap, tmesh$n)$x)
  grad <- gradient_of_linpoint(linpoint, smesh, tmesh)#
  prior.precision <- initial_Q
  
  cgen <- define.cgeneric.loggrow.model(linpoint, smesh, tmesh, step.size,
                                        prior.mean, prior.precision,
                                        priors = NULL, grad = grad,
                                        initial.growth = growth, initial.carry.cap = log(carry.cap), 
                                        initial.move.const = movement, initial.log.sigma = log(sigma), 
                                        debug = NULL)
  
  #mu_mat <- INLAtools::cgeneric_mu(cgen, theta)  #doesn't work, only gives dimension??
  
  a.func <- function(growth,carry.cap, linpoint){
    #print("Calcualting a")
    return(growth*exp(linpoint)/carry.cap)
  }
  r.vector <- function(growth,carry.cap,move.const,linpoint,grad){
    mag.grad.sq <- rowSums(grad*grad) #magnitude squared
    return(growth*exp(linpoint)*(linpoint-1)/carry.cap+ growth - move.const*mag.grad.sq )
  }
  fT <- function(a_array,movement, CinvG){
    return(movement*CinvG + Matrix::Diagonal(smesh$n, 1/step.size + a_array))
  }
  interpret.theta = function() {
    return(list(growth = exp(theta[1L]),
                carry.cap = exp(theta[2L]),
                move.const = exp(theta[3L]), 
                sigma = exp(theta[4L])))
  }
  mu = function(){
    out <- Matrix::Matrix(NA, nrow = smesh$n*tmesh$n, ncol = 1)
    par <- interpret.theta()
    r <- c(prior.mean, r.vector(par$growth, par$carry.cap, par$move.const, linpoint, grad)[-(1:smesh$n)])
    a_full <- a.func(par$growth, par$carry.cap, linpoint)
    out[1:smesh$n, 1] <- prior.mean
    fem.matrice <- fm_fem(smesh)
    CinvG <- Matrix::solve(fem.matrice$c1, fem.matrice$g1)
    for(t in 1:timesteps){
      fmat <- fT(a_full[t*smesh$n + 1:smesh$n], par$move.const, CinvG)
      out[t*smesh$n + 1:smesh$n,1] <- Matrix::solve(fmat, r[t*smesh$n + 1:smesh$n] + out[(t-1)*smesh$n + 1:smesh$n,1])
    }
    return(out)
  }
  Q_mat <- INLAtools::cgeneric_Q(cgen, theta)
  mu_mat = mu()
  #generate field
  if(debug) print("generating field")
  field <- data.frame(field = inla.qsample(1, Q_mat, mu = mu_mat)[, 1])
  field$time <- rep(0:timesteps, each = smesh$n)
  expand_for_plot <- function(i){
    animal_tempsf <- expand.grid(
      easting = seq(boundaries[1],boundaries[2], by = 0.01),
      northing = seq(boundaries[1],boundaries[2], by = 0.01))
    animal_tempsf <- dplyr::mutate(sf::st_as_sf(animal_tempsf, coords = c("easting", "northing")),
                                   time = i)
    animal_tempsf$field <- fmesher::fm_evaluate(
      smesh,
      loc = animal_tempsf,
      field = field$field[field$time == i])
    return(animal_tempsf)
  }
  expanded <- parallel::mclapply(0:timesteps, expand_for_plot,  mc.cores = ncores)
  animal <- do.call(rbind, expanded)
  bnd_inner <- sf::st_as_sf(inlabru::spoly(data.frame(easting = c(boundaries[1],boundaries[2],boundaries[2],boundaries[1]), 
                                                      northing = c(boundaries[1], boundaries[1], boundaries[2], boundaries[2]))))
  if(debug) print("Sampling")
  if(sample.type == "Normal"){
    if(is.null(obs.sd) | is.null(npoints)){
      warning("obs.sd and npoints must be defined")
    }
    points.to.sample <- sample(unique(sf::st_filter(animal,bnd_inner)$geometry),
                               npoints)
    animal_obs <- filter(animal, geometry %in% points.to.sample) %>% 
      dplyr::mutate(obs = rnorm(npoints*(tmesh$n), field, obs.sd))
  } else if(sample.type == "Bernoulli"){
    points.to.sample <- sample(unique(sf::st_filter(animal,bnd_inner)$geometry),
                               npoints)
    animal_obs <- filter(animal, geometry %in% points.to.sample) %>% 
      dplyr::mutate(obs = rbinom(npoints*(tmesh$n), 1, plogis(obs.prob + field)), 
                    survey = rep(1, npoints*tmesh$n))
    if(nsurv>1){
      for(i in 2:nsurv){
        animal_obs <- rbind(animal_obs, 
                            filter(animal, geometry %in% points.to.sample) %>% 
                              dplyr::mutate(obs = rbinom(npoints*(tmesh$n), 1, plogis(obs.prob+field)), 
                                            survey = rep(i, npoints*tmesh$n)))
      }
    }
  }else if(sample.type == "LGCP"){
    points <- st_as_sf(sample.lgcp(smesh, field$field[field$time == 0], samplers = bnd_inner),
                       coords = c("x","y"))
    animal_obs <- dplyr::mutate(points, time = 0)
    
    for(i in 1:timesteps){
      pointsi <- dplyr::mutate(st_as_sf(sample.lgcp(smesh, field$field[field$time == i], samplers = bnd_inner),
                                        coords = c("x","y")),
                               time = i)
      animal_obs <- rbind(animal_obs, pointsi)
      rm(pointsi)
    }
    
  }else if(sample.type == "Poisson"){
    if(is.null(npoints)){
      warning("obs.sd and npoints must be defined")
    }
    points.to.sample <- sample(unique(sf::st_filter(animal,bnd_inner)$geometry),
                               npoints)
    animal_obs <- filter(animal, geometry %in% points.to.sample) %>% 
      dplyr::mutate(obs = rpois(npoints*(tmesh$n), exp(field)))
  }else{
    print("Sampling type not recognised")
    animal_obs = 0
  }
  return(list(animal = animal[animal$time != 0,],field = field[field$time !=0,],
              animal_obs = animal_obs[animal_obs$time != 0,], mesh = smesh))
}
