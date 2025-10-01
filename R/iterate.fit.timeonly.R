#'Fit a non spatial logistic growth model with Poisson observations
#' @description
#' Iteratively fits the spatial logistic growth model. 
#' @param data spatial dataframe of gaussian observations, with column y = log(observed value)
#' @param tmesh mesh over equally spaced time points created with \code{fm_mesh_1d}
#' @param prior.mean estimated mean for first year of data.
#' @param prior.precision uncertainty for estimated mean of first year of data. 
#' @param max.iter maximum iterations to attempt
#' @param gamma dampening parameter for update rule
#' @param stop.crit stopping criteria for linearisation point update rule. Stop updating if mean(abs(new_linearisation_point-old_linearisation_point))<=stop.crit
#' @param priors named list of prior parameters, named \code{cc} (carrying capacity), \code{growth}, \code{sigma}. 
#' Each is a two element vector containing the mean and variance for parameter. 
#' @param initial.linpoint Optional. Starting guess for the linearisation point. If NULL, will be estimated within function
#' @param initial.growth,initial.carry.cap,initial.log.sigma Starting values for the growth, 
#' \emph{log} carrying capacity, movement constant and \emph{log} standard deviation
#' @param verbose logical supplied to INLA
#' @returns list containing final model fit, number of iterations \code{n}, matrix of all past linearisation points and list of all past model fits.  
#'@export
iterate.timeonly <- function(data, tmesh,nsurvey, step.size, prior.mean,
                                     prior.precision, max.iter = 100,gamma = 0.75,stop.crit = 0.05,
                                     priors = NULL,initial.linpoint = NULL, initial.growth=1, 
                                     initial.carry.cap=0.05, initial.log.sigma = log(1.5), 
                                     verbose = F){
  #browser()
  if(is.null(initial.linpoint)){
    initial.linpoint <- log(logit.nest(exp(prior.mean), initial.growth, exp(initial.carry.cap), tmesh$n)$x)
  }
  if(!is.matrix(initial.linpoint)) initial.linpoint <- as.matrix(initial.linpoint, ncol = 1)
  fit.list <- list()
  #Set up initial model
  log_growth_model <- define.loggrow.time.model(linpoint = initial.linpoint, tmesh = tmesh, step.size = step.size, 
                                                        prior.mean = prior.mean,
                                                        prior.precision = prior.precision, priors = priors,
                                                        initial.growth = initial.growth, 
                                                        initial.carry.cap = initial.carry.cap,
                                                        initial.log.sigma = initial.log.sigma)
  fit <- bru(y ~ loggrow(time, 
                         model = log_growth_model)-1,
             data = data, domain = list(time = tmesh),
             family = "poisson", options = list(verbose = verbose))
  fit.list[[1]]<-fit
  print("First fitting finished")
  n.nodes <- fit$misc$configs$nconfig
  nodes <- data.frame(log.prob=rep(NA,n.nodes))
  mat_list <- list()
  mean_list <- list()
  for(i in 1:n.nodes){
    nodes[i,]<- c(fit$misc$configs$config[[i]]$log.posterior)
    mat_list[[i]] <- fit$misc$configs$config[[i]]$Q[1:(tmesh$n), 1:(tmesh$n)]
    mean_list[[i]] <- fit$misc$configs$config[[i]]$improved.mean[1:(tmesh$n)]
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
  #print(new.mean)
  new.linpoint <- (1-gamma)*initial.linpoint +gamma*new.mean
  #Check that this linpoint isn't so extreme that it will cause issues
  #plot(new.linpoint)
  lp.mat <- cbind(initial.linpoint,new.linpoint)
  n <- 2
  #print(fit$summary.hyperpar$mean)
  
  #Iterate the updates
  while(n < max.iter & mean(abs(lp.mat[,n]-lp.mat[,n-1]))>stop.crit){
    log_growth_model <- define.loggrow.time.model(linpoint = as.vector(new.linpoint),
                                                          tmesh = tmesh, step.size = step.size, 
                                                          prior.mean = prior.mean,
                                                          prior.precision = prior.precision, priors = priors,
                                                          #initial.growth = fit$summary.hyperpar$mean[1], 
                                                          #initial.carry.cap = fit$summary.hyperpar$mean[2],
                                                          #initial.log.sigma = fit$summary.hyperpar$mean[3])
                                                          initial.growth = initial.growth, 
                                                          initial.carry.cap = initial.carry.cap,
                                                          initial.log.sigma = initial.log.sigma)
    print("Defined new model")
    fit <- bru(y ~ loggrow(time, 
                           model = log_growth_model)-1,
               data = data, domain = list(time = tmesh),
               family = "poisson", options = list(verbose = verbose))
    print(paste("Fitted new model", n))
    n.nodes <- fit$misc$configs$nconfig
    if(!is.numeric(n.nodes)){
      print("Failed to fit, trying again")
      fit <- bru(y ~ loggrow(time, 
                             model = log_growth_model)-1,
                 data = data, domain = list(time = tmesh),
                 family = "poisson", options = list(verbose = verbose))
      if(!is.numeric(fit$misc$configs$nconfig)){
        print("Failed again, returning model output")
        return(list(new.linpoint = new.linpoint,fit = fit, past.linpoints = lp.mat, fit.list = fit.list))
      }
      n.nodes <- fit$misc$configs$nconfig
    }
    fit.list[[n]]<-fit
    nodes <- data.frame(log.prob=rep(NA,n.nodes))
    mat_list <- list()
    mean_list <- list()
    for(i in 1:n.nodes){
      nodes[i,]<- c(fit$misc$configs$config[[i]]$log.posterior)
      mat_list[[i]] <- fit$misc$configs$config[[i]]$Q[1:(tmesh$n),1:(tmesh$n)]
      mean_list[[i]] <- fit$misc$configs$config[[i]]$improved.mean[1:(tmesh$n)]
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
    new.linpoint <- (1-gamma)*lp.mat[,n] +gamma*new.mean
    #plot(new.linpoint, main = paste("Linearisation point", n))
    lp.mat <- cbind(lp.mat,new.linpoint)
    print("Updated linpoint")
    n <- n+1
  }
  log_growth_model <- define.loggrow.time.model(linpoint = as.vector(new.linpoint),
                                                        tmesh = tmesh, step.size = step.size, 
                                                        prior.mean = prior.mean,
                                                        prior.precision = prior.precision, priors = priors,
                                                        #initial.growth = fit$summary.hyperpar$mean[1], 
                                                        #initial.carry.cap = fit$summary.hyperpar$mean[2],
                                                        #initial.log.sigma = fit$summary.hyperpar$mean[3])
                                                        initial.growth = initial.growth, 
                                                        initial.carry.cap = initial.carry.cap,
                                                        initial.log.sigma = initial.log.sigma)
  print("Fitting final model")
  final.fit <- bru(y ~ loggrow(time, 
                               model = log_growth_model)-1,
                   data = data, domain = list(time = tmesh),
                   family = "poisson", options = list(verbose = verbose))
  return(list(fit = final.fit, n = n, linpoints = lp.mat, fit.list = fit.list))
}