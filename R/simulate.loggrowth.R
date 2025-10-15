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
                               max.edge = 0.05){
  browser()
  #functions needed
  a.func <- function(growth,carry.cap, linpoint){
    #print("Calcualting a")
    return(growth*exp(linpoint)/carry.cap)
  }
  #growth, carry.cap, move.const = theta params to be est
  #step.size = difference in time between lin points, known
  #linpoint = list of linearisation point vectors 
  #smesh = space mesh built with fmesher, tmesh = time mesh
  L.matrix <- function(growth,carry.cap,move.const,step.size, linpoint, smesh, tmesh){
    #print("Calcualting Lmat")
    ns <- smesh$n
    nt <- tmesh$n
    a<- a.func(growth,carry.cap, linpoint)
    a[1:ns] <- 1
    a.mat <- Matrix::Diagonal(ns*nt,a)
    subdiag <- Matrix::bandSparse(nt*ns, k = -ns, diagonals = list(rep(-1/step.size, (nt - 1)*ns)))
    fem.matrices <- fmesher::fm_fem(smesh)
    CinvG <- Matrix::solve(fem.matrices$c1, fem.matrices$g1)
    main.diag <- Matrix::kronecker(Matrix::Diagonal(nt, c(0,rep(1, nt-1))), 
                                   Matrix::Diagonal(ns, 1/(step.size))+ move.const*CinvG)
    #print(diag(main.diag + subdiag + a.mat))
    return(Matrix::drop0(main.diag + subdiag + a.mat, tol = 1e-100))
  }
  r.vector <- function(growth,carry.cap,move.const,linpoint,grad){
    mag.grad.sq <- rowSums(grad*grad) #magnitude squared
    return(growth*exp(linpoint)*(linpoint-1)/carry.cap+ growth - move.const*mag.grad.sq )
  }
  interpret.theta = function() {
    return(list(growth = theta[1L],
                carry.cap = exp(theta[2L]),
                move.const = theta[3L], 
                sigma = exp(theta[4L])))
  }
  Q = function(){
    #browser()
    #print("Calcualting Q")
    par = interpret.theta()
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, par$move.const,step.size, linpoint, smesh, tmesh)
    noiseonly = Matrix::Diagonal(smesh$n*(tmesh$n-1), 1/(par$sigma*step.size)**2)
    noise.precision = Matrix::bdiag(list(prior.precision, noiseonly))
    output = Matrix::crossprod(Lmat, noise.precision %*% Lmat)
    #print(output[smesh$n:(smesh$n +10),smesh$n:(smesh$n +10)])
    return(Matrix::drop0(output, 1e-100))
  }
  mu = function(){
    #browser()
    #print("Calcualting mu")
    #if(class(theta)!="numeric"){
    #  theta <- initial()
    #}
    par = interpret.theta()
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, par$move.const, step.size, linpoint, smesh, tmesh)
    r = c(prior.mean, r.vector(par$growth, par$carry.cap, par$move.const, linpoint, grad)[-(1:smesh$n)])
    Lmat_det <- Matrix::det(Lmat)
    if(!is.nan(Lmat_det)) {
      if(abs(Lmat_det) <= .Machine$double.eps){ #if close to singular use
        #print(det(crossprod(Lmat,Lmat)))
        mu = Matrix::solve(crossprod(Lmat,Lmat),crossprod(Lmat,r)) #more stable form of solve(lmat,r)
        mu= as.vector(mu)
        #print("Trick version")
      }else{
        mu = Matrix::solve(Lmat,r)
        #print("Default Solve")
      }}else{
        print("There's some NaNs going on?")
        mu = NA
      }
    #print(mean(mu))
    return(mu)
  }
  
  #set up for simulation
  corners <- c(boundaries[1] - 2*initial.range, boundaries[2]+2*initial.range)
  bnd_extended <- inlabru::spoly(data.frame(easting = c(corners[1], corners[2],corners[2],corners[1]), 
                                            northing = c(corners[1], corners[1],corners[2],corners[2])))
  smesh <- fmesher::fm_mesh_2d_inla(boundary = bnd_extended, max.edge = max.edge)
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
    #summary(exp(prior.mean))
    print("Defining model")
  }
  #components needed for model
  initial.growth <- growth
  initial.carry.cap <- log(carry.cap)
  initial.move.const <- movement
  initial.log.sigma <- log(sigma)
  theta <- c(initial.growth, initial.carry.cap, initial.move.const, initial.log.sigma)
  linpoint <- log(logit.nest(exp(prior.mean), growth, exp(carry.cap), tmesh$n)$x)
  #plot(exp(linpoint))
  grad <- gradient_of_linpoint(linpoint, smesh, tmesh)#
  prior.precision <- initial_Q
  if(debug) print("Calculating precision")
  Q_mat <-Q()
  if(debug) print("Calculating mean")
  mu_mat <- mu()
  #if(debug) plot(exp(mu_mat))
  
  #generate field
  if(debug) print("generating field")
  field <- data.frame(field = inla.qsample(1, Q_mat, mu = mu_mat)[, 1])
  field$time <- rep(0:timesteps, each = smesh$n)
  expand_for_plot <- function(i){
    animal_tempsf <- expand.grid(
      easting = seq(corners[1],corners[2], by = 0.01),
      northing = seq(corners[1],corners[2], by = 0.01))
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
    points.to.sample <- sample(unique(sf::st_filter(animal,bnd_inner)$geometry),
                               npoints)
    animal_obs <- filter(animal, geometry %in% points.to.sample) %>% 
      dplyr::mutate(obs = rnorm(npoints*(tmesh$n), field, obs.sd))
  }
  if(sample.type == "LGCP"){
    simulate_obs <- function(i){
      samp_animal <- sample.lgcp(smesh, 
                                 loglambda = field$field[field$time == i],
                                 samplers = bnd_inner)
      samp_animal <- sf::st_as_sf(samp_animal, coords = c("x","y"))
      samp_animal_df <- dplyr::mutate(samp_animal, time = i)
      return(samp_animal_df)
    }
    observations <- parallel::mclapply(0:timesteps, simulate_obs,  mc.cores =  ncores)
    animal_obs <- do.call(rbind, observations)
    #remove edge effects
    #animal_obs <- st_as_sf(animal_obs, coords = c("x","y"))
  }
  if(sample.type == "Bernoulli"){
    points.to.sample <- sample(unique(sf::st_filter(animal,bnd_inner)$geometry),
                               npoints)
    animal_obs <- filter(animal, geometry %in% points.to.sample) %>% 
      dplyr::mutate(obs = rbinom(npoints*(tmesh$n), 1, plogis(obs.prob*exp(field))))
  }
  return(list(animal = animal[animal$time !=0,],field = field[field$time != 0,],
              animal_obs = animal_obs[animal_obs$time != 0,]))
}
