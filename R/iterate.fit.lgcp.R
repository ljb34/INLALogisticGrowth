#'Fit a logistic growth model with LGCP observations
#' @description
#' Iteratively fits the spatial logistic growth model. 
#' @param data spatial dataframe  observations
#' @param smesh spatial mesh created with \code{fm_mesh_2d_inla}
#' @param tmesh mesh over equally spaced time points created with \code{fm_mesh_1d}
#' @param samplers sampling area for supplying to inlabru
#' @param prior.mean estimated mean for first year of data. Helper function for calculating coming soon
#' @param prior.variance uncertainty for estimated mean of first year of data. Helper function for calculating coming soon
#' @param max.iter maximum iterations to attempt
#' @param gamma dampening parameter for update rule
#' @param stop.crit stopping criteria for linearisation point update rule. Stop updating if mean(abs(new_linearisation_point-old_linearisation_point))<=stop.crit
#' @param priors named list of prior parameters, named \code{cc} (carrying capacity), \code{growth}, \code{move}, \code{sigma}. 
#' Each is a two element vector containing the mean and variance for parameter. 
#' @param initial.linpoint Optional. Starting guess for the linearisation point. If NULL, will be estimated within function
#' @param initial.growth,initial.carry.cap,initial.move.const,initial.log.sigma Starting values for the growth, 
#' \emph{log} carrying capacity, movement constant and \emph{log} standard deviation
#' @param verbose logical supplied to INLA
#' @returns list containing final model fit, number of iterations \code{n}, matrix of all past linearisation points and list of all past model fits.  
#'@export
iterate.fit.lgcp <- function(data, smesh, tmesh, samplers,prior.mean,
                             prior.variance, max.iter = 100,gamma = 0.5,stop.crit = 0.05,
                             priors = NULL, initial.linpoint = NULL, initial.growth=1, 
                             initial.carry.cap=100, initial.move.const = 1, initial.log.sigma = log(1.5),
                             verbose = F){
  #browser()
  step.size = (tmesh$interval[2]-tmesh$interval[1])/(tmesh$n-1) #calculate step size. -1 in denom due to fence post problem 
  if(is.null(initial.linpoint)){
    initial.linpoint <- log(logit.nest(exp(prior.mean), exp(initial.growth), initial.carry.cap, tmesh$n)$x)
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
             family = "cp", options = list(verbose = verbose, control.inla = list(control.vb=list(emergency=30))))
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
  nodes <- dplyr::mutate(nodes, weight = exp(log.prob)) %>%
    dplyr::mutate(weight.prob = weight/sum(weight))
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
                                             prior.precision = prior.precision, priors = priors,
                                             initial.growth = initial.growth, 
                                             initial.carry.cap = initial.carry.cap,
                                             initial.move.const = initial.move.const,
                                             initial.log.sigma = initial.log.sigma)

    print("Defined new model")
    fit <- bru(geometry + time ~ loggrow(list(space = geometry, time = time), 
                                         model = log_growth_model, 
                                         n = smesh$n*tmesh$n) -1,
               data = data, domain = list(geometry = smesh,time = tmesh),
               samplers = samplers,
               family = "cp", options = list(verbose = verbose, control.inla = list(control.vb=list(emergency=30))))
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
                 family = "cp", options = list(verbose = verbose, control.inla = list(control.vb=list(emergency=300))))
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
    nodes <- dplyr::mutate(nodes, weight = exp(log.prob)) %>%
      dplyr::mutate(weight.prob = weight/sum(weight))
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
                                           prior.precision = prior.precision, priors = priors,
                                           initial.growth = initial.growth, 
                                           initial.carry.cap = initial.carry.cap,
                                           initial.move.const = initial.move.const,
                                           initial.log.sigma = initial.log.sigma)
  print("Defined final model")
  final.fit <- bru(geometry + time ~ loggrow(list(space = geometry, time = time), 
                                             model = log_growth_model, 
                                             n = smesh$n*tmesh$n) -1,
                   data = data, domain = list(geometry = smesh,time = tmesh),
                   samplers = samplers,
                   family = "cp", options = list(verbose = verbose))
  return(list(fit = final.fit, n = n, linpoints = lp.mat, fit_list = fit.list))
}