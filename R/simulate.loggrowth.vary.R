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
  
  Q = function(par){
    #browser()
    #print("Calcualting Q")
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, par$move.const,step.size, linpoint, smesh, tmesh)
    mats <- fmesher::fm_fem(smesh)
    g <- mean(par$move.const)
    noiseblock <- (mats$c1 + g*mats$g1)*(par$sigma**2)/step.size
    noiseonly = Matrix::bdiag(replicate(tmesh$n-1, noiseblock, simplify = FALSE))
    noise.precision = Matrix::bdiag(list(prior.precision, noiseonly))
    output = Matrix::crossprod(Lmat, noise.precision %*% Lmat)
    #print(output[smesh$n:(smesh$n +10),smesh$n:(smesh$n +10)])
    return(Matrix::drop0(output, 1e-100))
  }
  mu = function(par){
    #browser()
    #print("Calcualting mu")
    #if(class(theta)!="numeric"){
    #  theta <- initial()
    #}
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, par$move.const, step.size, linpoint, smesh, tmesh)
    r = c(prior.mean, r.vector(par$growth, par$carry.cap, par$move.const, linpoint, grad)[-(1:smesh$n)])
    Lmat_det <- Matrix::det(Lmat)
    if(!is.nan(Lmat_det)) {
      if(abs(Lmat_det) <= .Machine$double.eps){ #if close to singular use
        #print(det(crossprod(Lmat,Lmat)))
        mu = Matrix::solve(crossprod(Lmat,Lmat),crossprod(Lmat,r)) #more stable form of solve(lmat,r)
        mu= as.vector(mu)
        print("Trick version")
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
  if(debug) print("Calculating precision")
  Q_mat <-Q(par)
  if(debug) print("Calculating mean")
  mu_mat <- mu(par)
  
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
              animal_obs = animal_obs[animal_obs$time != 0,], mesh = smesh, covaraiates = cov.grid))
}

test <- simulate_loggrowth_vary(log(1.5), 0.1, log(1000),0.15, 0.2, 0.05, 20,
                                    0.2, 1,
                                    550,0.25, 0.05,4, debug = T, max.edge = 0.25,
                                same.cov = F)
ggplot()+gg(test$animal, aes(fill = exp(field)), geom = "tile")+
  facet_wrap(~time)+
  scale_fill_viridis_c()

ggplot()+gg(test$covaraiates, aes(fill = movement), geom = "tile")

