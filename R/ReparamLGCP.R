library(dplyr)
library(INLA)
library(fmesher)
library(inlabru)
library(sf)
library(sp)
library(Matrix)

# RGeneric set up ---------------------------------------------------------
log.growth.rgeneric =  function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL){ 
  envir = parent.env(environment()) #gets extra parameters (linpoint etc.) from definition data
  library(Matrix)
  library(fmesher)
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
    a.mat <- Diagonal(ns*nt,a)
    
    subdiag <- kronecker(bandSparse(nt, k = -1, diagonals = list(rep(1, nt - 1))),
                         Diagonal(ns, -1/(step.size)))
    fem.matrices <- fm_fem(smesh)
    CinvG <- solve(fem.matrices$c1, fem.matrices$g1)
    main.diag <- kronecker(Diagonal(nt, c(0,rep(1, nt-1))), 
                           Diagonal(ns, 1/(step.size))+ move.const*CinvG)
    #print(diag(main.diag + subdiag + a.mat))
    return(main.diag + subdiag + a.mat)
  }
  r.vector <- function(growth,carry.cap,move.const,linpoint,smesh, tmesh){
    #find 2 nearest neighbours to approximate gradient
    #browser()
    #print("Calcualting rvector")
    ns = smesh$n; nt = tmesh$n
    coords <- smesh$loc[,c(1,2)]
    distances <- as.matrix(dist(coords, upper = T))
    near.neighbours <- apply(distances, 2, order)[2:4,]
    grad <- matrix(nrow = ns*nt, ncol = 2)
    for(i in 1:ns){
      diffmat <- matrix(c(coords[near.neighbours[1,i],1]- coords[i,1], 
                          coords[near.neighbours[1,i],2]- coords[i,2],
                          coords[near.neighbours[2,i],1]- coords[i,1], 
                          coords[near.neighbours[2,i],2]- coords[i,2]),
                        byrow = T, nrow = 2)
      diffmat[which(abs(diffmat) < .Machine$double.eps, arr.ind = T)] <- 0
      for(t in 0:(nt-1)){
        if(abs(det(diffmat))<=.Machine$double.eps){ # if both nearest neighbours are exactly horizontal or both vertical from point, then go to 
          #1st and 3rd near neighbours
          diffmat2 <- matrix(c(coords[near.neighbours[1,i],1]- coords[i,1], 
                               coords[near.neighbours[1,i],2]- coords[i,2],
                               coords[near.neighbours[3,i],1]- coords[i,1], 
                               coords[near.neighbours[3,i],2]- coords[i,2]),
                             byrow = T, nrow = 2)
          diffmat[which(abs(diffmat) < .Machine$double.eps, arr.ind = T)] <- 0
          grad[t*ns+i,] <- solve(diffmat2,
                                 c(linpoint[near.neighbours[1,i]+t*ns] - linpoint[i + t*ns],
                                   linpoint[near.neighbours[2,i]+t*ns]- linpoint[i + t*ns]))
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
  interpret.theta = function() {
    return(list(growth = theta[1L],
                carry.cap = exp(theta[2L]),
                move.const = theta[3L], 
                sigma = exp(theta[4L])))
  }
  
  graph = function() {
    return (Q())
  }
  Q = function(){
    #print("Calcualting Q")
    par = interpret.theta()
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, par$move.const,step.size, linpoint, smesh, tmesh)
    noiseonly = Diagonal(smesh$n*(tmesh$n-1), (par$sigma*step.size)**2)
    noise.variance = bdiag(list(prior.variance, noiseonly))
    output = crossprod(Lmat, solve(noise.variance, Lmat))
    #print(output[smesh$n:(smesh$n +10),smesh$n:(smesh$n +10)])
    return(output)
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
    r = c(prior.mean, r.vector(par$growth, par$carry.cap, par$move.const, linpoint, smesh, tmesh)[-(1:smesh$n)])
    #print(det(Lmat))
    if(!is.nan(det(Lmat))) {
      if(abs(det(Lmat)) <= .Machine$double.eps|(is.infinite(det(Lmat)) & !is.infinite(det(crossprod(Lmat,Lmat))))){ #if close to singular use
        #print(det(crossprod(Lmat,Lmat)))
        mu = solve(crossprod(Lmat,Lmat),crossprod(Lmat,r)) #more stable form of solve(lmat,r)
        mu= as.vector(mu)
        #print("Trick version")
      }else{
        mu = solve(Lmat,r)
        #print("Default Solve")
      }}else{
        print("There's some NaNs going on?")
        mu = NA
      }
    #print(mean(mu))
    return(mu)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function(){#can change params to make user specified
    #print("Calcualting logprior")
    par = interpret.theta()
    #print(par)
    if(!is.null(priors)) warning("Parameters missing for priors")
    val = dnorm(par$carry.cap, mean = priors$cc[1], sd = priors$cc[2], log = T)+
      dnorm(par$growth, mean = priors$growth[1], sd = priors$growth[2], log = T)+
      dnorm(par$move.const,mean = priors$move[1], sd = priors$move[2], log = T)+ 
      dnorm(par$sigma, mean = priors$sigma[1], sd = priors$sigma[2], log = T)
    return(val)
  }
  initial = function(){#can change params to make user specified
    if(is.null(initial.growth)) initial.growth = 1
    if(is.null(initial.carry.cap)) initial.carry.cap = 100
    if(is.null(initial.move.const)) initial.move.const = 1
    if(is.null(initial.log.sigma)) initial.log.sigma = log(5)
    return(c(initial.growth, initial.carry.cap, initial.move.const, initial.log.sigma))
  }
  quit = function() {
    return(invisible())
  }
  if (is.null(theta)) theta = initial()
  if (length(theta) == 0) theta = initial()
  val = do.call(match.arg(cmd), args = list())
  return(val)
}

define.loggrow.model <- function(linpoint, smesh, tmesh, step.size,
                                 prior.mean, prior.variance, priors = NULL,
                                 initial.growth = NULL, initial.carry.cap = NULL, 
                                 initial.move.const = NULL, initial.log.sigma = NULL){
  
  the_model <- inla.rgeneric.define(log.growth.rgeneric, 
                                    linpoint = linpoint, 
                                    smesh = smesh,tmesh = tmesh, step.size = step.size, 
                                    prior.mean = prior.mean, priors = priors,
                                    prior.variance = prior.variance, 
                                    initial.growth = initial.growth, 
                                    initial.carry.cap = initial.carry.cap, 
                                    initial.move.const = initial.move.const, 
                                    initial.log.sigma = initial.log.sigma)
  class(the_model) <- c("log_growth_model", class(the_model))
  the_model[["smesh"]] <- smesh
  the_model[["tmesh"]] <- tmesh
  return(the_model)
}


bru_get_mapper.log_growth_model <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  inlabru::bru_mapper_multi(list(
    space = inlabru::bru_mapper(model[["smesh"]]),
    time = inlabru::bru_mapper(model[["tmesh"]], indexed = TRUE)
  ))
}


# Iterated ----------------------------------------------------------------
iterate.fit.lgcp <- function(data, smesh, tmesh, samplers,prior.mean,
                             prior.variance, max.iter = 100,gamma = 0.5,stop.crit = 0.05,
                             priors = NULL, initial.linpoint = NULL, initial.growth=1, 
                             initial.carry.cap=100, initial.move.const = 1, initial.log.sigma = log(1.5),
                             verbose = F){
  #browser()
  step.size = (tmesh$interval[2]-tmesh$interval[1])/(tmesh$n-1) #calculate step size. -1 in denom due to fence post problem 
  if(is.null(initial.linpoint)){
    initial.linpoint <- log(logit.nest(exp(prior.mean), initial.growth, exp(initial.carry.cap), tmesh$n)$x)
  }
  if(!is.matrix(initial.linpoint)) initial.linpoint <- as.matrix(initial.linpoint, ncol = 1)
  #Set up initial model
  fit.list <- list()
  log_growth_model <- define.loggrow.model(linpoint = initial.linpoint, 
                                           smesh = smesh,tmesh = tmesh, step.size = step.size, 
                                           prior.mean = prior.mean,
                                           prior.variance = prior.variance, priors = priors,
                                           initial.growth = initial.growth, 
                                           initial.carry.cap = initial.carry.cap,
                                           initial.move.const = initial.move.const,
                                           initial.log.sigma = initial.log.sigma)
  fit <- bru(geometry + time ~ loggrow(list(space = geometry, time = time), 
                                       model = log_growth_model, 
                                       n = smesh$n*tmesh$n) -1,
             data = data, domain = list(geometry = smesh,time = tmesh),
             samplers = samplers,
             family = "cp", options = list(verbose = verbose))
  fit.list[[1]]<-fit
  print("First fitting finished")
  n.nodes <- fit$misc$configs$nconfig
  nodes <- data.frame(log.prob=rep(NA,n.nodes))
  #mat_list <- list()
  mean_list <- list()
  for(i in 1:n.nodes){
    nodes[i,]<- fit$misc$configs$config[[i]]$log.posterior
    #mat_list[[i]] <- fit$misc$configs$config[[i]]$Q[1:(smesh$n*tmesh$n), 1:(smesh$n*tmesh$n)]
    mean_list[[i]] <- fit$misc$configs$config[[i]]$improved.mean[1:(smesh$n*tmesh$n)]
  }
  nodes <- mutate(nodes, weight = exp(log.prob)) %>%
    mutate(weight.prob = weight/sum(weight))
  #Old rule- in theory faster but gives some extreme changes
  #P <- Reduce("+", Map(function(m, w) m * w, mat_list, nodes$weight.prob))
  #weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
  #b <- Reduce("+", Map(function(m,w) m%*%w, mat_list,weighted.means))
  #new.linpoint <- (1-gamma)*initial.linpoint +gamma*solve(P,b)
  
  #New update rule
  weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
  new.mean <- Reduce("+", weighted.means)
  new.linpoint <- (1-gamma)*initial.linpoint +gamma*new.mean
  print("Calcualted new linpoint")
  lp.mat <- cbind(initial.linpoint,new.linpoint)
  n <- 2
  #print(fit$summary.hyperpar$mean)
  while(n < max.iter & mean(abs(lp.mat[,n]-lp.mat[,n-1]))>stop.crit){
    log_growth_model <- define.loggrow.model(linpoint = as.vector(new.linpoint), 
                                             smesh = smesh,tmesh = tmesh, step.size = step.size, 
                                             prior.mean = prior.mean,
                                             prior.variance = prior.variance, priors = priors,
                                             initial.growth = fit$summary.hyperpar$mean[1], 
                                             initial.carry.cap = fit$summary.hyperpar$mean[2],
                                             initial.move.const = fit$summary.hyperpar$mean[3],
                                             initial.log.sigma = fit$summary.hyperpar$mean[4])
    print("Defined new model")
    fit <- bru(geometry + time ~ loggrow(list(space = geometry, time = time), 
                                         model = log_growth_model, 
                                         n = smesh$n*tmesh$n) -1,
               data = data, domain = list(geometry = smesh,time = tmesh),
               samplers = samplers,
               family = "cp", options = list(verbose = verbose))
    print(paste("Fitted new model", n))
    fit.list[[n]]<-fit
    n.nodes <- fit$misc$configs$nconfig
    if(!is.numeric(n.nodes)){
      print("Failed to fit, trying again")
      fit <- bru(geometry + time ~ loggrow(list(space = geometry, time = time), 
                                           model = log_growth_model, 
                                           n = smesh$n*tmesh$n) -1,
                 data = data, domain = list(geometry = smesh,time = tmesh),
                 samplers = samplers,
                 family = "cp", options = list(verbose = verbose))
      n.nodes <- fit$misc$configs$nconfig
      if(!is.numeric(fit$misc$configs$nconfig)){
        print("Failed again, returning model output")
        return(list(new.linpoint = new.linpoint,fit = fit, past.linpoints = lp.mat, fit.list = fit.list))
      }
      fit.list[[n]]<-fit
    }
    nodes <- data.frame(log.prob=rep(NA,n.nodes))
    #mat_list <- list()
    mean_list <- list()
    for(i in 1:n.nodes){
      nodes[i,]<- fit$misc$configs$config[[i]]$log.posterior
      #mat_list[[i]] <- fit$misc$configs$config[[i]]$Q[1:(smesh$n*tmesh$n),1:(smesh$n*tmesh$n)]
      mean_list[[i]] <- fit$misc$configs$config[[i]]$improved.mean[1:(smesh$n*tmesh$n)]
    }
    nodes <- mutate(nodes, weight = exp(log.prob)) %>%
      mutate(weight.prob = weight/sum(weight))
    #P <- Reduce("+", Map(function(m, w) m * w, mat_list, nodes$weight.prob))
    #weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
    #b <- Reduce("+", Map(function(m,w) m%*%w, mat_list,weighted.means))
    #new.linpoint <- (1-gamma)*lp.mat[,n] +gamma*solve(P,b)
    
    #New update rule
    weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
    new.mean <- Reduce("+", weighted.means)
    new.linpoint <- (1-gamma)*lp.mat[,n-1] +gamma*new.mean
    
    lp.mat <- cbind(lp.mat,new.linpoint)
    print("Updated linpoint")
    n <- n+1
  }
  log_growth_model <- define.loggrow.model(linpoint = as.vector(new.linpoint), 
                                           smesh = smesh,tmesh = tmesh, step.size = step.size, 
                                           prior.mean = prior.mean,
                                           prior.variance = prior.variance, priors = priors,
                                           initial.growth = fit$summary.hyperpar$mean[1], 
                                           initial.carry.cap = fit$summary.hyperpar$mean[2],
                                           initial.move.const = fit$summary.hyperpar$mean[3],
                                           initial.log.sigma = fit$summary.hyperpar$mean[4])
  print("Defined final model")
  final.fit <- bru(geometry + time ~ loggrow(list(space = geometry, time = time), 
                                             model = log_growth_model, 
                                             n = smesh$n*tmesh$n) -1,
                   data = data, domain = list(geometry = smesh,time = tmesh),
                   samplers = samplers,
                   family = "cp", options = list(verbose = verbose))
  return(list(fit = final.fit, n = n, linpoints = lp.mat, fit_list = fit.list))
}


# Simulate data -----------------------------------------------------------
logit.growth <- function(x,r,k){
  x[x<0]<-0
  xnew <- x*exp(r*(1-(x/k)))
  return(xnew)
}
logit.nest <- function(x0,r,k,n){
  df <- data.frame(x = x0, time = rep(1, length(x0)))
  for(i in 2:n){
    df <- rbind(df,
                data.frame(x = logit.growth(df$x[df$time == i-1], r,k), 
                           time = rep(i, length(x0))))
  }
  return(df)
}
simulate.loggrowth<- function(growth, k, movement, sigma, 
                              initial, timesteps, npoints = NULL, obs.sd=NULL, 
                              sample.type = "LGCP", ncores = 1,
                              boundaries = c(0,1)){
  #browser()
  #Set up boundary and mesh
  corners <- c(boundaries[1] - movement, boundaries[2]+movement)
  bnd_extended <- spoly(data.frame(easting = c(corners[1], corners[2],corners[2],corners[1]), 
                                   northing = c(corners[1], corners[1],corners[2],corners[2])))
  mesh_extended <- fm_mesh_2d_inla(boundary = bnd_extended, max.edge = 0.05)
  mesh_time <- fm_mesh_1d(loc = 1:timesteps)
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
    animal_tempsf <- mutate(sf::st_as_sf(animal_tempsf, coords = c("easting", "northing")),
                            time = i)
    animal_tempsf$field <- fm_evaluate(
      mesh_extended,
      loc = animal_tempsf,
      field = animal_field$field[animal_field$time == i])
    return(animal_tempsf)
  }
  expanded <- parallel::mclapply(1:timesteps, expand_for_plot,  mc.cores = ncores)
  animal <- do.call(rbind, expanded)
  bnd_inner <- st_as_sf(spoly(data.frame(easting = c(boundaries[1],boundaries[2],boundaries[2],boundaries[1]), 
                                         northing = c(boundaries[1], boundaries[1], boundaries[2], boundaries[2]))))
  if(sample.type == "Normal"){
    points.to.sample <- sample(unique(st_filter(animal,bnd_inner)$geometry),
                               npoints)
    animal_obs <- filter(animal, geometry %in% points.to.sample) %>% 
      mutate(obs = rnorm(npoints*(timesteps+1), field, obs.sd))
  }
  if(sample.type == "LGCP"){
    animal_field$field[animal_field$field <0] <- 0.00001
    simulate_obs <- function(i){
      samp_animal <- sample.lgcp(mesh_extended, 
                                 loglambda = log(animal_field$field[animal_field$time == i]),
                                 samplers = bnd_inner)
      samp_animal <- st_as_sf(samp_animal, coords = c("x","y"))
      samp_animal_df <- mutate(samp_animal, time = i)
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

# The good stuff ----------------------------------------------------------
out.lgcp <- simulate.loggrowth(growth = 1, k = 150, movement = 1, sigma = 20,
                               initial = 50,timesteps = 4,sample.type = "LGCP")

#fit initial year
bnd <- spoly(data.frame(easting = c(0,1,1,0), northing = c(0,0,1,1)))
mesh_obs <- fm_mesh_2d_inla(boundary = bnd,
                            max.edge = c(0.2,1))
bnd <- st_as_sf(bnd)
mesh_time <- fm_mesh_1d(loc = 1:4)
matern <- inla.spde2.pcmatern(mesh_obs,
                              prior.sigma = c(0.1, 0.1),
                              prior.range = c(0.1, 0.1))
cmp <- geometry ~ smooth(geometry, model = matern) + 
  initial(1,model = "linear", mean.linear = length(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,]))-1

fit0 <- bru(cmp, out.lgcp$animal_obs[out.lgcp$animal_obs$time == 1,],domain = list(geometry = mesh_obs),
            family = "cp",samplers = bnd)
#Find fitted values on mesh points
library(stringr)
index <- min(which(str_sub(rownames(fit0$summary.fitted.values),8,8)!= "A"))

initial.variance <- Diagonal(mesh_obs$n, (fit0$summary.fixed$sd**2)+(fit0$summary.fitted.values$sd[index-1 +1:mesh_obs$n]**2))
#fit other years
priors <- list(cc = c(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,]),50),
               growth = c(1,1),move = c(1,1),sigma = c(log(20),1))
iterated.fit.lgcp <- iterate.fit.lgcp(data = out.lgcp$animal_obs, smesh = mesh_obs, tmesh = mesh_time,
                                      samplers = bnd,prior.mean = fit0$summary.fixed$mean +fit0$summary.fitted.values$mean[index-1 +1:mesh_obs$n],
                                      prior.variance = initial.variance, priors = priors,
                                      max.iter = 100,gamma = 0.5,
                                      stop.crit = 0.01,
                                      initial.linpoint = NULL, initial.growth = 0.8, 
                                      initial.carry.cap = log(nrow(out.lgcp$animal_obs[out.lgcp$animal_obs$time == 4,])),
                                      verbose = T)
iterated.fit.lgcp$data <- out.lgcp
saveRDS(iterated.fit.lgcp, "LogGrowth/LGCPfit_reparam.RData")
