#'Functions for simulating growth
simulate.loggrowth2<- function(growth, k, movement, sigma, 
                               initial, timesteps, npoints = NULL, obs.sd=NULL, 
                               sample.type = "LGCP", ncores = 1,
                               boundaries = c(0,1)){
  browser()
  corners <- c(boundaries[1] - movement, boundaries[2]+movement)
  bnd_extended <- spoly(data.frame(easting = c(corners[1], corners[2],corners[2],corners[1]), 
                                   northing = c(corners[1], corners[1],corners[2],corners[2])))
  mesh_space <- fm_mesh_2d_inla(boundary = bnd_extended, max.edge = 0.3)
  mesh_time <- fm_mesh_1d(loc = 1:timesteps)
  
  matrices <- inla.logrowth.mean.precision(mesh_space, mesh_time, c(growth, log(k), movement, log(sigma)),
                                       initial)
  
  field <- inla.qsample(1, matrices$Q, mu = matrices$mean)[, 1]
  field$time <- rep(1:timesteps, each = mesh_space$time)
  expand_for_plot <- function(i){
    animal_tempsf <- expand.grid(
      easting = seq(corners[1],corners[2], by = 0.01),
      northing = seq(corners[1],corners[2], by = 0.01))
    animal_tempsf <- mutate(sf::st_as_sf(animal_tempsf, coords = c("easting", "northing")),
                            time = i)
    animal_tempsf$field <- fm_evaluate(
      mesh_extended,
      loc = animal_tempsf,
      field = field$field[field$time == i])
    return(animal_tempsf)
  }
  expanded <- parallel::mclapply(1:timesteps, expand_for_plot,  mc.cores = ncores)
  animal <- do.call(rbind, expanded)
  
  return(list(field = field, animal = animal))
}
#'Get precision matrix

inla.logrowth.mean.precision <- function(smesh, tmesh, theta, initial){
  browser()
  movement <- theta[3]
  sigma <- exp(theta[4])
  linpoint <- logit.nest(rep(initial,smesh$n), theta[1], exp(theta[2]), tmesh$n)$x
  matern <-
    inla.spde2.pcmatern(smesh,
                        prior.sigma = c(0.1, 0.1),
                        prior.range = c(0.1, 0.1))
  initial_Q <- inla.spde.precision(matern,
                                   theta = log(c(movement, sigma)))
  logit_model <- define.loggrow.model(linpoint = linpoint, smesh = smesh, tmesh = tmesh,
                                      step.size = 1, prior.mean = rep(initial, smesh$n),
                                      prior.variance = solve(initial_Q))
  precision <- inla.rgeneric.q(cmd = "Q", rmodel = logit_model, theta = theta)
  mu <- inla.rgeneric.q(cmd = "mu", rmodel = logit_model, theta = theta)
  return(list(Q = precision, mean = mu))
}





bnd_extended <- spoly(data.frame(easting = c(0,2,2,0), 
                                               northing = c(0,0,2,2)))
mesh_extended <- fm_mesh_2d_inla(boundary = bnd_extended, max.edge = 0.25)
mesh_time <- fm_mesh_1d(loc = 1:4)
test.p<- inla.logrowth.precision(mesh_extended, mesh_time, c(1,log(1/20),1,1), 10)

test.out <- simulate.loggrowth2(0.75,500,1,1,200,4)

