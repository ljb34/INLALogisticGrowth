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

simulate_loggrowth<- function(growth, k, movement, sigma, 
                              initial, timesteps, npoints = NULL, obs.sd=NULL, 
                              sample.type = "LGCP", ncores = 1,
                              boundaries = c(0,1)){
  #browser()
  #Set up boundary and mesh
  corners <- c(boundaries[1] - movement, boundaries[2]+movement)
  bnd_extended <- inlabru::spoly(data.frame(easting = c(corners[1], corners[2],corners[2],corners[1]), 
                                   northing = c(corners[1], corners[1],corners[2],corners[2])))
  mesh_extended <- fmesher::fm_mesh_2d_inla(boundary = bnd_extended, max.edge = 0.05)
  mesh_time <- fmesher::fm_mesh_1d(loc = 1:timesteps)
  #animal initial field
  matern_extended <-
    inla.spde2.pcmatern(mesh_extended,
                        prior.sigma = c(0.1, 0.1),
                        prior.range = c(0.1, 0.1))
  initial_Q <- inla.spde.precision(matern_extended,
                                   theta = log(c(movement, sigma)))
  initial_field <- inla.qsample(1, initial_Q, mu = rep(initial, nrow(initial_Q)))[, 1]
  animal_field <- data.frame(time = 1, field = initial_field)
  for(i in 2:timesteps){
    mu <- logit.growth(animal_field$field[animal_field$time == i-1],
                       growth,k)
    animal_field<-rbind(animal_field, 
                        data.frame(time = i, 
                                   field = inla.qsample(1, initial_Q, mu = mu)[, 1]))
  }
  #Expand for plots
  expand_for_plot <- function(i){
    animal_tempsf <- expand.grid(
      easting = seq(corners[1],corners[2], by = 0.01),
      northing = seq(corners[1],corners[2], by = 0.01))
    animal_tempsf <- dplyr::mutate(sf::st_as_sf(animal_tempsf, coords = c("easting", "northing")),
                            time = i)
    animal_tempsf$field <- fmesher::fm_evaluate(
      mesh_extended,
      loc = animal_tempsf,
      field = animal_field$field[animal_field$time == i])
    return(animal_tempsf)
  }
  expanded <- parallel::mclapply(1:timesteps, expand_for_plot,  mc.cores = ncores)
  animal <- do.call(rbind, expanded)
  bnd_inner <- sf::st_as_sf(inlabru::spoly(data.frame(easting = c(boundaries[1],boundaries[2],boundaries[2],boundaries[1]), 
                                         northing = c(boundaries[1], boundaries[1], boundaries[2], boundaries[2]))))
  if(sample.type == "Normal"){
    points.to.sample <- sample(unique(sf::st_filter(animal,bnd_inner)$geometry),
                               npoints)
    animal_obs <- filter(animal, geometry %in% points.to.sample) %>% 
      dplyr::mutate(obs = rnorm(npoints*(timesteps), field, obs.sd))
  }
  if(sample.type == "LGCP"){
    animal_field$field[animal_field$field <0] <- 0.00001
    simulate_obs <- function(i){
      samp_animal <- sample.lgcp(mesh_extended, 
                                 loglambda = log(animal_field$field[animal_field$time == i]),
                                 samplers = bnd_inner)
      samp_animal <- sf::st_as_sf(samp_animal, coords = c("x","y"))
      samp_animal_df <- dplyr::mutate(samp_animal, time = i)
      return(samp_animal_df)
    }
    observations <- parallel::mclapply(1:timesteps, simulate_obs,  mc.cores =  ncores)
    animal_obs <- do.call(rbind, observations)
    #remove edge effects
    #animal_obs <- st_as_sf(animal_obs, coords = c("x","y"))
  }
  return(list(animal = animal,animal_field = animal_field,
              animal_obs = animal_obs))
}

