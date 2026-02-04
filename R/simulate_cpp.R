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
#'field defined on mesh nodes for debugging purposes (field). 
#'
#'@export
simulate_cpp <- function(growth, carry.cap, movement, sigma, 
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
  hex_points <- fm_hexagon_lattice(bnd = bnd_extended, edge_len = 0.9*max.edge[1])
  smesh <- fmesher::fm_mesh_2d_inla(loc = hex_points, boundary = bnd_extended,
                                    max.edge = max.edge,
                                    offset = c(-0.01, (boundaries[2]-boundaries[1])/2))
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
  linpoint <- log(logit.nest(exp(prior.mean), growth, carry.cap, tmesh$n)$x)
  grad <- gradient_of_linpoint(linpoint, smesh, tmesh)#
  prior.precision <- as(initial_Q, "dgCMatrix")
  fem_matrices <- fmesher::fm_fem(smesh)
  CinvG <-Matrix::solve(fem_matrices$c1, fem_matrices$g1)
  fem_matrices$c1 <- as(fem_matrices$c1, "dgCMatrix")
  fem_matrices$g1 <- as(fem_matrices$g1, "dgCMatrix")
  CinvG <- as(CinvG,"dgCMatrix")
  
  Q_mat <- Q_sparse_cpp(growth, carry.cap, movement, sigma, timestep = 1, linpoint, 
                    smesh$n, tmesh$n, fem_matrices$c1, fem_matrices$g1, 
                    CinvG, prior.precision)
  mu_mat <- mu_sparse_cpp(growth, carry.cap, movement, timestep = 1, linpoint, 
                      grad, prior.mean, smesh$n, tmesh$n, CinvG)
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
    lambda_func <- function(x,y){
      points_sf <- st_as_sf(data.frame(x = x, y = y), coords = c("x","y"))
      field_values <- fmesher::fm_evaluate(
        smesh,
        loc = points_sf,
        field = field$field[field$time == 0])
      return(exp(field_values))
    }
    spatstat_sim <- spatstat.random::rpoispp(lambda = lambda_func, win = spatstat.geom::owin(boundaries,boundaries))
    spatstat_df <- data.frame(x = spatstat_sim$x, y = spatstat_sim$y, time = rep(0, length(spatstat_sim$x)))
    animal_obs <- st_as_sf(spatstat_df, coords = c("x","y"))
    
    for(i in 1:timesteps){
      lambda_func_i <- function(x,y){
        points_sf <- st_as_sf(data.frame(x = x, y = y), coords = c("x","y"))
        field_values <- fmesher::fm_evaluate(
          smesh,
          loc = points_sf,
          field = field$field[field$time == i])
        return(exp(field_values))
      }
      spatstat_sim_i <- spatstat.random::rpoispp(lambda = lambda_func_i, win = spatstat.geom::owin(boundaries,boundaries))
      spatstat_df_i <- data.frame(x = spatstat_sim_i$x, y = spatstat_sim_i$y, time = rep(i, length(spatstat_sim_i$x)))
      spatstat_sf_i <- st_as_sf(spatstat_df_i, coords = c("x","y"))
      animal_obs <- rbind(animal_obs, spatstat_sf_i)
      rm(spatstat_sim_i, spatstat_df_i, spatstat_sf_i)
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
  return(list(animal = animal[animal$time != 0,],field = field[field$time != 0,],
              animal_obs = animal_obs[animal_obs$time != 0,], mesh = smesh))
}

