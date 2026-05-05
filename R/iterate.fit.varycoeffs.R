
#' @param covariates = named list of covariates
#'@export 
iterate.fit.varycoeffs <- function(formula, data,family, smesh, tmesh, samplers,prior.mean,
                               prior.precision,growth.formula = ~1, 
                               carry.formula = ~1, move.formula = ~1, covariates,
                               max.iter = 100,gamma = 0.5,stop.crit = 0.05,
                               priors = NULL, initial.linpoint = NULL, initial.growth=0.5, 
                               initial.carry.cap=1000, initial.move.const = 0.5, initial.log.sigma = log(1.5),
                               update.rule = 2, debug = F, options = NULL, saveall = T,
                               weights = NULL,domain = NULL, early.stop = F, ...){
  #browser()
  step.size = (tmesh$interval[2]-tmesh$interval[1])/(tmesh$n-1) #calculate step size. -1 in denom due to fence post problem 
  if(is.null(initial.linpoint)){
    initial.linpoint <- log(logit.nest(exp(prior.mean), exp(initial.growth[1]), exp(initial.carry.cap[1]), tmesh$n)$x)
  }
  if(!is.matrix(initial.linpoint)) initial.linpoint <- as.matrix(initial.linpoint, ncol = 1)
  fit_list <- list()
  if(is.null(domain)){
    domain = list(geometry = smesh, time = tmesh)
  }
  #Covariates
  vars_growth <- all.vars(growth.formula)
  vars_carry  <- all.vars(carry.formula)
  vars_move   <- all.vars(move.formula)
  
  all_vars <- unique(c(vars_growth, vars_carry, vars_move))
  ns <- smesh$n
  nt <- tmesh$n
  
  coords <- smesh$loc[, 1:2]
  time_index <- 1:nt
  
  mesh_df <- expand.grid(space = 1:ns, time = 1:nt)
  mesh_df$x <- coords[mesh_df$space, 1]
  mesh_df$y <- coords[mesh_df$space, 2]
  
  for (v in all_vars) {
    r <- covariates[[v]]
    
    if (inherits(r, "SpatRaster")) {
      
      if (terra::nlyr(r) == 1) {
        # spatial-only → repeat across time
        vals <- terra::extract(r, coords)[,2]
        mesh_df[[v]] <- vals[mesh_df$space]
        
      } else if (terra::nlyr(r) == nt) {
        # space-time raster (layers = time)
        
        vals_mat <- matrix(NA, ns, nt)
        
        for (t in 1:nt) {
          vals_mat[, t] <- terra::extract(r[[t]], coords)[,2]
        }
        
        mesh_df[[v]] <- vals_mat[
          cbind(mesh_df$space, mesh_df$time)
        ]
        
      } else {
        stop(paste("Raster", v, "has incompatible number of layers"))
      }
      
    } else {
      stop(paste("Covariate", v, "is not a SpatRaster"))
    }
  }
  growth_cov <- as.vector(model.matrix(growth.formula, data = mesh_df))
  carry_cov  <- as.vector(model.matrix(carry.formula,  data = mesh_df))
  move_cov   <- as.vector(model.matrix(move.formula,   data = mesh_df))
  #Set up initial model
  log_growth_model <- define.varying.cgeneric.loggrow.model(linpoint = initial.linpoint, 
                                                      smesh = smesh,tmesh = tmesh, step.size = step.size, 
                                                      prior.mean = prior.mean,
                                                      prior.precision = prior.precision,
                                                      growth.formula = growth.formula,
                                                      carry.formula = carry.formula,
                                                      move.formula = move.formula,
                                                      growth_cov = growth_cov,
                                                      carry_cov = carry_cov,
                                                      move_cov = move_cov,
                                                      priors = priors,
                                                      initial.growth = initial.growth, 
                                                      initial.carry.cap = initial.carry.cap,
                                                      initial.move.const = initial.move.const,
                                                      initial.log.sigma = initial.log.sigma, debug = debug)
  
  new.cmp <- update(formula, . ~ . + loggrow(list(space = geometry, time = time),
                                             model = log_growth_model, n = smesh$n * tmesh$n))
  environment(new.cmp) <- environment()
  
  fit <- bru(new.cmp,
             data = data, domain = domain,
             samplers = samplers,
             family = family, options = options,
             weights = weights, ...)
  print("Fitted new model 1")
  if(saveall){
    fit_list[[1]]<-fit
  }else{
    fit_list <- fit
  }
  if(early.stop){
    return(list(fit = fit, linpoints = initial.linpoint))
  }
  n.nodes <- fit$misc$configs$nconfig
  nodes <- data.frame(log.prob=rep(NA,n.nodes))
  mat_list <- list()
  mean_list <- list()
  for(i in 1:n.nodes){
    nodes[i,]<- fit$misc$configs$config[[i]]$log.posterior
    Q <- fit$misc$configs$config[[i]]$Q[1:(smesh$n*tmesh$n),1:(smesh$n*tmesh$n)]
    dQ <- Matrix::diag(Q)
    Q <- Q + Matrix::t(Q)
    Matrix::diag(Q) <- dQ
    mat_list[[i]] <- Q
    mean_list[[i]] <- fit$misc$configs$config[[i]]$improved.mean[1:(smesh$n*tmesh$n)]
  }
  nodes <- dplyr::mutate(nodes, weight = exp(log.prob)) %>%
    dplyr::mutate(weight.prob = weight/sum(weight))
  if(update.rule == 2){
    print("Performing type II update")
    #Type II update
    P <- Reduce("+", Map(function(m, w) m * w, mat_list, nodes$weight.prob))
    weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
    b <- Reduce("+", Map(function(m,w) m%*%w, mat_list,weighted.means))
    new.linpoint <- (1-gamma)*initial.linpoint +gamma*Matrix::solve(P,b)
  }else{
    #Type I
    weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
    new.mean <- Reduce("+", weighted.means)
    print(length(new.mean))
    print(length(initial.linpoint))
    new.linpoint <- (1-gamma)*initial.linpoint +gamma*new.mean
  }
  
  print("Calculated new linpoint")
  lp.mat <- cbind(initial.linpoint,new.linpoint)
  n <- 2
  #print(fit$summary.hyperpar$mean)
  while(n < max.iter & mean(abs(lp.mat[,n]-lp.mat[,n-1]))>stop.crit){
    
    log_growth_model <- define.varying.cgeneric.loggrow.model(linpoint = as.vector(new.linpoint), 
                                                        smesh = smesh,tmesh = tmesh, step.size = step.size, 
                                                        prior.mean = prior.mean,
                                                        prior.precision = prior.precision,
                                                        growth.formula = growth.formula,
                                                        carry.formula = carry.formula,
                                                        move.formula = move.formula,
                                                        growth_cov = growth_cov,
                                                        carry_cov = carry_cov,
                                                        move_cov = move_cov,
                                                        priors = priors,
                                                        initial.growth = initial.growth, 
                                                        initial.carry.cap = initial.carry.cap,
                                                        initial.move.const = initial.move.const,
                                                        initial.log.sigma = initial.log.sigma, debug = debug)
    
    
    print("Defined new model")
    new.cmp <- update(formula, . ~ . + loggrow(list(space = geometry, time = time),
                                               model = log_growth_model, n = smesh$n * tmesh$n))
    environment(new.cmp) <- environment()
    fit <- bru(new.cmp,
               data = data, domain = domain,
               samplers = samplers,
               family = family, options = options,
               weights = weights, ...)
    print(paste("Fitted new model", n))
    if(saveall){
      fit_list[[n]]<-fit
    }else{
      fit_list <- fit
    }
    n.nodes <- fit$misc$configs$nconfig
    if(!is.numeric(n.nodes)){
      print("Failed to fit, trying again")
      fit <- bru(geometry + time ~ loggrow(list(space = geometry, time = time), 
                                           model = log_growth_model, 
                                           n = smesh$n*tmesh$n) -1,
                 data = data, domain = domain,
                 samplers = samplers,
                 family = "cp", options = options,
                 weights = weights, ... )
      n.nodes <- fit$misc$configs$nconfig
      if(!is.numeric(fit$misc$configs$nconfig)){
        print("Failed again, returning model output")
        if(saveall){
          fit_list[[n]]<-fit
        } else{
          fit_list <- fit
        }
        return(list(new.linpoint = new.linpoint,fit = fit, linpoints = lp.mat, fit_list = fit_list))
      }
    }
    nodes <- data.frame(log.prob=rep(NA,n.nodes))
    mat_list <- list()
    mean_list <- list()
    for(i in 1:n.nodes){
      nodes[i,]<- fit$misc$configs$config[[i]]$log.posterior
      Q <- fit$misc$configs$config[[i]]$Q[1:(smesh$n*tmesh$n),1:(smesh$n*tmesh$n)]
      dQ <- Matrix::diag(Q)
      Q <- Q + Matrix::t(Q)
      Matrix::diag(Q) <- dQ
      mat_list[[i]] <- Q
      mean_list[[i]] <- fit$misc$configs$config[[i]]$improved.mean[1:(smesh$n*tmesh$n)]
    }
    nodes <- dplyr::mutate(nodes, weight = exp(log.prob)) %>%
      dplyr::mutate(weight.prob = weight/sum(weight))
    if(update.rule == 2){
      #Type II update
      P <- Reduce("+", Map(function(m, w) m * w, mat_list, nodes$weight.prob))
      weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
      b <- Reduce("+", Map(function(m,w) m%*%w, mat_list,weighted.means))
      new.linpoint <- (1-gamma)*initial.linpoint +gamma*Matrix::solve(P,b)
    }else{
      #Type I
      weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
      new.mean <- Reduce("+", weighted.means)
      new.linpoint <- (1-gamma)*initial.linpoint +gamma*new.mean
    }
    
    lp.mat <- cbind(lp.mat,new.linpoint)
    print("Updated linpoint")
    print(summary(exp(new.linpoint)))
    n <- n+1
    if(saveall){
      fit_list[[n]]<-fit
    } else{
      fit_list <- fit
    }
  }

  log_growth_model <- define.varying.cgeneric.loggrow.model(linpoint = as.vector(new.linpoint), 
                                                      smesh = smesh,tmesh = tmesh, step.size = step.size, 
                                                      prior.mean = prior.mean,
                                                      prior.precision = prior.precision,
                                                      growth.formula = growth.formula,
                                                      carry.formula = carry.formula,
                                                      move.formula = move.formula,
                                                      growth_cov = growth_cov,
                                                      carry_cov = carry_cov,
                                                      move_cov = move_cov,
                                                      priors = priors,
                                                      initial.growth = initial.growth, 
                                                      initial.carry.cap = initial.carry.cap,
                                                      initial.move.const = initial.move.const,
                                                      initial.log.sigma = initial.log.sigma, debug = debug)
  
  print("Defined final model")
  new.cmp <- update(formula, . ~ . + loggrow(list(space = geometry, time = time),
                                             model = log_growth_model, n = smesh$n * tmesh$n))
  environment(new.cmp) <- environment()
  final.fit <- bru(new.cmp,
                   data = data, domain = domain,
                   samplers = samplers,
                   family = family, options = options,
                   weights = weights, ...)
  return(list(fit = final.fit, n = n, linpoints = lp.mat, fit_list = fit_list))
}
