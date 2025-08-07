#'Functions for simulating growth
#'@export
simulate.loggrowth2<- function(growth, k, movement, sigma, 
                               initial, timesteps, npoints = NULL, obs.sd=NULL, 
                               sample.type = "LGCP", ncores = 1,
                               boundaries = c(0,1)){
  #browser()
  corners <- c(boundaries[1] - movement, boundaries[2]+movement)
  bnd_extended <- inlabru::spoly(data.frame(easting = c(corners[1], corners[2],corners[2],corners[1]), 
                                   northing = c(corners[1], corners[1],corners[2],corners[2])))
  mesh_space <- fmesher::fm_mesh_2d_inla(boundary = bnd_extended, max.edge = 0.2)
  mesh_time <- fmesher::fm_mesh_1d(loc = 1:timesteps)
  print("set up finished, calculating mean and precision")
  matrices <- inla.logrowth.mean.precision(mesh_space, mesh_time, c(growth, log(k), movement, log(sigma)),
                                       initial)
  print("Calculated mean and precision, generating field")
  field <- data.frame(field = inla.qsample(1, matrices$Q, mu = matrices$mean)[, 1])
  field$time <- rep(1:timesteps, each = mesh_space$n)
  expand_for_plot <- function(i){
    animal_tempsf <- expand.grid(
      easting = seq(corners[1],corners[2], by = 0.01),
      northing = seq(corners[1],corners[2], by = 0.01))
    animal_tempsf <- dplyr::mutate(sf::st_as_sf(animal_tempsf, coords = c("easting", "northing")),
                            time = i)
    animal_tempsf$field <- fmesher::fm_evaluate(
      mesh_space,
      loc = animal_tempsf,
      field = field$field[field$time == i])
    return(animal_tempsf)
  }
  expanded <- parallel::mclapply(1:timesteps, expand_for_plot,  mc.cores = ncores)
  animal <- do.call(rbind, expanded)
  bnd_inner <- sf::st_as_sf(inlabru::spoly(data.frame(easting = c(boundaries[1],boundaries[2],boundaries[2],boundaries[1]), 
                                                      northing = c(boundaries[1], boundaries[1], boundaries[2], boundaries[2]))))
  print("Sampling")
  if(sample.type == "Normal"){
    points.to.sample <- sample(unique(sf::st_filter(animal,bnd_inner)$geometry),
                               npoints)
    animal_obs <- filter(animal, geometry %in% points.to.sample) %>% 
      mutate(obs = rnorm(npoints*(timesteps), field, obs.sd))
  }
  if(sample.type == "LGCP"){
    field$field[field$field <0] <- 0.00001
    simulate_obs <- function(i){
      samp_animal <- sample.lgcp(mesh_space, 
                                 loglambda = field$field[field$time == i],
                                 samplers = bnd_inner)
      samp_animal <- sf::st_as_sf(samp_animal, coords = c("x","y"))
      samp_animal_df <- mutate(samp_animal, time = i)
      return(samp_animal_df)
    }
    observations <- parallel::mclapply(1:timesteps, simulate_obs,  mc.cores =  ncores)
    animal_obs <- do.call(rbind, observations)
    #remove edge effects
    #animal_obs <- st_as_sf(animal_obs, coords = c("x","y"))
  }
  return(list(animal = animal,field = field,
              animal_obs = animal_obs))
}
#'Get precision matrix
#'@export
inla.logrowth.mean.precision <- function(smesh, tmesh, theta, initial){
  #browser()
  growth <- theta[1]
  carry.cap <- exp(theta[2])
  movement <- theta[3]
  sigma <- exp(theta[4])
  step.size = (tmesh$interval[2]-tmesh$interval[1])/(tmesh$n-1)
  matern <-
    inla.spde2.pcmatern(smesh,
                        prior.sigma = c(0.1, 0.1),
                        prior.range = c(0.1, 0.1))
  initial_Q <- inla.spde.precision(matern,
                                   theta = c(movement, log(1/sigma)))
  initial.mean <- inla.qsample(1, initial_Q, mu = rep(log(initial), nrow(initial_Q)))[,1]
  linpoint <- log(logit.nest(exp(initial.mean), growth, carry.cap, tmesh$n)$x)
  logit_model <- define.loggrow.model(linpoint = linpoint, smesh = smesh, tmesh = tmesh,
                                      step.size = step.size, prior.mean = initial.mean,
                                      prior.variance = solve(initial_Q))
  precision <- inla.rgeneric.q(cmd = "Q", rmodel = logit_model, theta = theta)
  mu <- mu(growth, carry.cap, movement, step.size, linpoint, smesh, tmesh, initial.mean)
  return(list(Q = precision, mean = mu))
}

#'L matrix
#'@export
L.matrix <- function(growth,carry.cap,move.const,step.size, linpoint, smesh, tmesh){
  print("Calcualting Lmat")
  #browser()
  ns <- smesh$n
  nt <- tmesh$n
  a<- a.func(growth,carry.cap, linpoint)
  a[1:ns] <- 1
  a.mat <- Matrix::Diagonal(ns*nt,a)
  
  subdiag <- Matrix::kronecker(Matrix::bandSparse(nt, k = -1, diagonals = list(rep(1, nt - 1))),
                               Matrix::Diagonal(ns, -1/(step.size)))
  fem.matrices <- fmesher::fm_fem(smesh)
  CinvG <- Matrix::solve(fem.matrices$c1, fem.matrices$g1)
  main.diag <- Matrix::kronecker(Matrix::Diagonal(nt, c(0,rep(1, nt-1))), 
                                 Matrix::Diagonal(ns, 1/(step.size))+ move.const*CinvG)
  #print(diag(main.diag + subdiag + a.mat))
  return(main.diag + subdiag + a.mat)
}
#'r vector
#'@export
r.vector <- function(growth,carry.cap,move.const,linpoint,smesh, tmesh){
  #find 2 nearest neighbours to approximate gradient
  #browser()
  print("Calcualting rvector")
  ns = smesh$n; nt = tmesh$n
  coords <- smesh$loc[,c(1,2)]
  distances <- as.matrix(dist(coords, upper = T))
  near.neighbours <- apply(distances, 2, order)[2:10,]
  grad <- matrix(nrow = ns*nt, ncol = 2)
  for(i in 1:ns){
    diffmat <- matrix(c(coords[near.neighbours[1,i],1]- coords[i,1], 
                        coords[near.neighbours[1,i],2]- coords[i,2],
                        coords[near.neighbours[2,i],1]- coords[i,1], 
                        coords[near.neighbours[2,i],2]- coords[i,2]),
                      byrow = T, nrow = 2)
    diffmat[which(abs(diffmat) < .Machine$double.eps, arr.ind = T)] <- 0
    for(t in 0:(nt-1)){
      if(abs(Matrix::det(diffmat))<=.Machine$double.eps){ # if both nearest neighbours are exactly horizontal or both vertical from point, then go to 
        #1st and 3rd near neighbours
        j <-3
        while(min(abs(c(coords[near.neighbours[1,i],1]-coords[near.neighbours[j,i],1],
                       coords[near.neighbours[1,i],2]-coords[near.neighbours[j,i],2]))) < .Machine$double.eps |
              coords[near.neighbours[1,i],1]%%coords[near.neighbours[j,i],1] == coords[near.neighbours[1,i],2]%%coords[near.neighbours[j,i],2]){
          j <- j+1
          if(j == 10){
            warning(paste("Mesh behaving strangely. All nearest points to point", i, "are exactly horizontal or vertical."))
          }
        }
        diffmat2 <- matrix(c(coords[near.neighbours[1,i],1]- coords[i,1], 
                             coords[near.neighbours[1,i],2]- coords[i,2],
                             coords[near.neighbours[j,i],1]- coords[i,1], 
                             coords[near.neighbours[j,i],2]- coords[i,2]),
                           byrow = T, nrow = 2)
        diffmat2[which(abs(diffmat2) < .Machine$double.eps, arr.ind = T)] <- 0
        grad[t*ns+i,] <- solve(diffmat2,
                               c(linpoint[near.neighbours[1,i]+t*ns] - linpoint[i + t*ns],
                                 linpoint[near.neighbours[j,i]+t*ns]- linpoint[i + t*ns]))
      }else{                  
        grad[t*ns+i,] <- solve(diffmat,
                               c(linpoint[near.neighbours[1,i]+t*ns] - linpoint[i + t*ns],
                                 linpoint[near.neighbours[2,i]+t*ns]- linpoint[i + t*ns]))
      }
    }
  }
  mag.grad.sq <- rowSums(grad*grad) #magnitude squared
  return(growth*exp(linpoint)*(linpoint-1)/carry.cap+ growth - move.const*mag.grad.sq )
}
#'mean
#'@export
mu = function(growth, carry.cap, move.const, step.size, linpoint, smesh, tmesh, prior.mean){
  #browser()
  print("Calculating mean")
  Lmat = L.matrix(growth, carry.cap, move.const, step.size, linpoint, smesh, tmesh)
  r = c(prior.mean, r.vector(growth, carry.cap, move.const, linpoint, smesh, tmesh)[-(1:smesh$n)])
  #print(det(Lmat))
  if(!is.nan(Matrix::det(Lmat))) {
    if(abs(Matrix::det(Lmat)) <= .Machine$double.eps|(is.infinite(Matrix::det(Lmat)) & !is.infinite(Matrix::det(crossprod(Lmat,Lmat))))){ #if close to singular use
      #print(det(crossprod(Lmat,Lmat)))
      mu = Matrix::solve(Matrix::crossprod(Lmat,Lmat),Matrix::crossprod(Lmat,r)) #more stable form of solve(lmat,r)
      mu= as.vector(mu)
      #print("Trick version")
    }else{
      mu = solve(Lmat,r)
      #print("Default Solve")
    }}else{
      warning("There's some NaNs going on?")
      mu = NA
    }
  #print(mean(mu))
  return(mu)
}
#' a helper function
#'@export
a.func <- function(growth,carry.cap, linpoint){
  #print("Calcualting a")
  return(growth*exp(linpoint)/carry.cap)
}

