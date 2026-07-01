#'@export
simulate_loggrowth_vary <- function(growth0, growth1, carry.cap0,carry.cap1, movement0, movement1, sigma,
                                    cov.range, cov.sigma,
                               initial.pop,initial.range, initial.sigma, 
                               timesteps, npoints = NULL, obs.sd=NULL,
                               obs.prob = NULL, same.cov = T,
                               sample.type = "LGCP", ncores = 1,
                               boundaries = c(0,1), debug = F,
                               max.edge = 0.05, nsurv = 3){
  #browser()
  #functions needed
  a.func <- function(growth,carry.cap, linpoint){
    #print("Calcualting a")
    return(growth*exp(linpoint)/carry.cap)
  }
  #growth, carry.cap, move.const = theta params to be est
  #step.size = difference in time between lin points, known
  #linpoint = list of linearisation point vectors 
  #smesh = space mesh built with fmesher, tmesh = time mesh
  r.vector <- function(growth,carry.cap,move.const,linpoint,grad){
    mag.grad.sq <- rowSums(grad*grad) #magnitude squared
    return(growth*exp(linpoint)*(linpoint-1)/carry.cap+ growth - move.const*mag.grad.sq )
  }
  
  fT <- function(a_array,movement, CinvG){
    return(movement*CinvG + Matrix::Diagonal(smesh$n, 1/step.size + a_array))
  }
  mu = function(){
    out <- Matrix::Matrix(NA, nrow = smesh$n*tmesh$n, ncol = 1)
    r <- c(prior.mean, r.vector(par$growth, par$carry.cap, par$move.const, linpoint, grad)[-(1:smesh$n)])
    a_full <- a.func(par$growth, par$carry.cap, linpoint)
    out[1:smesh$n, 1] <- prior.mean
    fem.matrice <- fm_fem(smesh)
    CinvG <- Matrix::solve(fem.matrice$c0, fem.matrice$g1)
    for(t in 1:timesteps){
      fmat <- fT(a_full[t*smesh$n + 1:smesh$n], par$move.const[t*smesh$n + 1:smesh$n], CinvG)
      out[t*smesh$n + 1:smesh$n,1] <- Matrix::solve(fmat, r[t*smesh$n + 1:smesh$n] + out[(t-1)*smesh$n + 1:smesh$n,1])
    }
    return(out)
  }
  
  Q = function(){
    #browser()
    out <- Matrix::Matrix(data = rep(0,(smesh$n*tmesh$n)**2), nrow = smesh$n*tmesh$n, ncol = smesh$n*tmesh$n)
    fem.matrice <- fm_fem(smesh)
    Qblock <- fem.matrice$c0 + par$move.const[1:smesh$n]*fem.matrice$g1
    CinvG <- Matrix::solve(fem.matrice$c0, fem.matrice$g1)
    a_full <- a.func(par$growth, par$carry.cap, linpoint)
    out[1:smesh$n, 1:smesh$n] <- initial_Q + (1/(par$sigma**3*step.size**2))*Qblock
    ft1 <- fT(a_full[smesh$n + 1:smesh$n], par$move.const[smesh$n + 1:smesh$n], CinvG)
    out[smesh$n + 1:smesh$n, 1:smesh$n] <- (-1/(par$sigma**2*step.size**2))*Qblock%*%ft1
    fmat <- ft1
    for(t in 1:(timesteps-1)){
      out[t*smesh$n + 1:smesh$n, (t-1)*smesh$n + 1:smesh$n] <- (-1/(par$sigma**2*step.size**2))*t(fmat)%*%Qblock
      out[t*smesh$n + 1:smesh$n, t*smesh$n + 1:smesh$n] <- (1/(par$sigma**2*step.size))*t(fmat)%*%Qblock%*%fmat
      ft1 <- fT(a_full[(t+1)*smesh$n + 1:smesh$n], par$move.const[(t+1)*smesh$n + 1:smesh$n], CinvG)
      out[t*smesh$n + 1:smesh$n, (t+1)*smesh$n + 1:smesh$n]<- (-1/(par$sigma**2*step.size**2))*ft1%*%Qblock
      fmat <- ft1
    }
    t = timesteps
    out[t*smesh$n + 1:smesh$n, (t-1)*smesh$n + 1:smesh$n] <- (-1/(par$sigma**2*step.size**2))*t(fmat)%*%Qblock
    out[t*smesh$n + 1:smesh$n, t*smesh$n + 1:smesh$n] <- (1/(par$sigma**2*step.size))*t(fmat)%*%Qblock%*%fmat
    return(out)
  }
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
  cov_mesh <- fmesher::fm_mesh_2d_inla( boundary = bnd_extended,
                                                max.edge = c(max.edge/3, max.edge),
                                                offset = c(-0.01, (boundaries[2]-boundaries[1])))
  cov_matern <-
    inla.spde2.pcmatern(cov_mesh,
                        prior.sigma = c(0.1, 0.1),
                        prior.range = c(0.1, 0.1))
  cov_Q <- inla.spde.precision(cov_matern,theta = log(c(cov.range, cov.sigma)))
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
    cov_mesh,
    loc = cov.grid,
    field = covariates$growth)
  cov.grid$carry.cap <- fmesher::fm_evaluate(
    cov_mesh,
    loc = cov.grid,
    field = covariates$carry.cap)
  cov.grid$movement <- fmesher::fm_evaluate(
    cov_mesh,
    loc = cov.grid,
    field = covariates$movement)
  
  mesh_pts <- sf::st_as_sf(
    data.frame(x = smesh$loc[,1], y = smesh$loc[,2]),
    coords = c("x", "y")
  )
  nn <- sf::st_nearest_feature(
    mesh_pts,
    cov.grid
  )
  
  mesh_pts$growth <- cov.grid$growth[nn]
  mesh_pts$carry.cap <- cov.grid$carry.cap[nn]
  mesh_pts$movement <- cov.grid$movement[nn]
  
  growth <- rep(exp(growth0 + growth1*mesh_pts$growth),timesteps+1)
  carry.cap <- rep(exp(carry.cap0 + carry.cap1*mesh_pts$carry.cap), timesteps+1)
  move.const <- rep(movement0 + movement1*mesh_pts$movement, timesteps+1)
  print(summary(carry.cap))
  print(summary(growth))
  par = list(growth = growth, carry.cap = carry.cap, move.const = move.const, sigma = sigma)
  print(length(prior.mean))
  print(length(growth)/(timesteps + 1))
  print(smesh$n)
  linpoint <- log(logit.nest(exp(prior.mean), growth[1:smesh$n], carry.cap[1:smesh$n], tmesh$n)$x)
  grad <- gradient_of_linpoint(linpoint, smesh, tmesh)#
  prior.precision <- initial_Q
  #browser()
  if(debug) print("Calculating precision")
  # cgen <- define.varying.cgeneric.loggrow.model(linpoint, smesh, tmesh, step.size,
  #                                       prior.mean, prior.precision, growth.formula = ~1 + growth,
  #                                       carry.formula = ~1 + carry.cap, move.formula = ~1 + movement,
  #                                       growth_cov = mesh_pts$growth, carry_cov = mesh_pts$carry.cap,
  #                                       move_cov = mesh_pts$movement,
  #                                       priors = NULL, grad = grad,
  #                                       initial.growth = growth, initial.carry.cap = log(carry.cap), 
  #                                       initial.move.const = move.const, initial.log.sigma = log(sigma), 
  #                                       debug = NULL)
  #Q_mat <- INLAtools::cgeneric_Q(cgen, theta = c(growth0, growth1, carry.cap0, carry.cap1, movement0, movement1, sigma))
  Q_mat <- Q()
  Qtest <- solve(Q_mat)
  if(debug) print("Calculating mean")
  mu_mat <- mu()
  summary(exp(mu_mat))
  #generate field
  if(debug) print("generating field")
  field <- data.frame(field = inla.qsample(1, Matrix::drop0(Q_mat), mu = mu_mat)[, 1])
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
  print(summary(animal))
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
  return(list(animal = animal[animal$time !=0,],field = field[field$time !=0,],
              animal_obs = animal_obs[animal_obs$time != 0,], mesh = smesh, covariates = cov.grid))
}



