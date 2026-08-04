iterate.timeonly.vary <- function(formula = ~-1, data, tmesh, step.size, prior.mean,
                             prior.precision, max.iter = 100,gamma = 0.75,stop.crit = 0.05,
                             priors = NULL, growth.formula = ~1, 
                             carry.formula = ~1, covariates = NULL,
                             initial.linpoint = NULL, initial.growth=1, 
                             initial.carry.cap=0.05, initial.sigma = log(1.5), 
                             verbose = F, family = "gaussian", domain = NULL){
  #browser()
  if(is.null(initial.linpoint)){
    initial.linpoint <- log(logit.nest(exp(prior.mean), initial.growth, exp(initial.carry.cap), tmesh$n)$x)
  }
  if(!is.matrix(initial.linpoint)) initial.linpoint <- as.matrix(initial.linpoint, ncol = 1)
  fit.list <- list()
  if(is.null(domain)){
    domain = list(time = seq(min(data$time), max(data$time), by = step.size))
  }
  #Covariates
  vars_growth <- all.vars(growth.formula)
  vars_carry  <- all.vars(carry.formula)
  if(attr(terms(growth.formula), "intercept")){
    growth_cov <- matrix(nrow = tmesh$n, ncol = length(vars_growth) + 1)
    growth_cov[,1] <- rep(1, tmesh$n)
    if(length(vars_growth) >=1){
      for(i in 1:length(vars_growth)){
        if(vars_growth[i] %in% names(covariates)){
          growth_cov[,1+i] <- covariates[[vars_growth[i]]]
        }else if(vars_growth[i] %in% names(data)){
          growth_cov[,1+i] <- data[, vars_growth[i]]
        } else{
          warning(paste("Couldn't find covariate in dataframe or covariates list ", vars_growth[i]))
        }
      }
    }
  } else if(length(vars_growth) >= 1){
    growth_cov <- matrix(nrow = tmesh$n, ncol = length(vars_growth))
    for(i in 1:length(vars_growth)){
      if(vars_growth[i] %in% names(covariates)){
        growth_cov[,i] <- covariates[[vars_growth[i]]]
      } else if(vars_growth[i] %in% names(data)){
        growth_cov[,i] <- data[, vars_growth[i]]
      } else {
        warning(paste("Couldn't find covariate in dataframe or covariates list ", vars_growth[i]))
      }
    }
  }
  if(attr(terms(carry.formula), "intercept")){
    carry_cov <- matrix(nrow = tmesh$n, ncol = length(vars_carry) + 1)
    carry_cov[,1] <- rep(1, tmesh$n)
    if(length(vars_carry) >=1){
      for(i in 1:length(vars_carry)){
         if(vars_carry[i] %in% names(covariates)){
          carry_cov[,1+i] <- covariates[[vars_carry[i]]]
        } else if(vars_carry[i] %in% names(data)){
          carry_cov[,1+i] <- data[, vars_carry[i]]
        }else {
          warning(paste("Couldn't find covariate in dataframe or covariates list ", vars_carry[i]))
        }
      }
    }
  } else if(length(vars_carry) >= 1){
    carry_cov <- matrix(nrow = tmesh$n, ncol = length(vars_carry))
    for(i in 1:length(vars_carry)){
      if(vars_carry[i] %in% names(data)){
        carry_cov[,i] <- data[, vars_carry[i]]
      } else if(vars_carry[i] %in% names(covariates)){
        carry_cov[,i] <- covariates[[vars_carry[i]]]
      } else{
        warning(paste("Couldn't find covariate in dataframe or covariates list ", vars_carry[i]))
      }
    }
  }
  #Set up initial model
  log_growth_model <- define.loggrow.vary.time.model(linpoint = initial.linpoint, tmesh = tmesh, step.size = step.size, 
                                                prior.mean = prior.mean,
                                                prior.precision = prior.precision, priors = priors,
                                                initial.growth = initial.growth, 
                                                initial.carry.cap = initial.carry.cap,
                                                initial.sigma = initial.sigma,
                                                growth.formula = growth.formula, 
                                                growth_cov = growth_cov,
                                                carry.formula = carry.formula,
                                                carry_cov = carry_cov)
  new.cmp <- update(formula, . ~ . + loggrow(time,
                                             model = log_growth_model))
  environment(new.cmp) <- environment()
  fit <- bru(new.cmp,
             data = data, domain = domain,
             family = family, options = list(verbose = verbose))
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
    log_growth_model <- define.loggrow.vary.time.model(linpoint = initial.linpoint, tmesh = tmesh, step.size = step.size, 
                                                       prior.mean = prior.mean,
                                                       prior.precision = prior.precision, priors = priors,
                                                       initial.growth = initial.growth, 
                                                       initial.carry.cap = initial.carry.cap,
                                                       initial.sigma = initial.sigma,
                                                       growth.formula = growth.formula, 
                                                       growth_cov = growth_cov,
                                                       carry.formula = carry.formula,
                                                       carry_cov = carry_cov)
    new.cmp <- update(formula, . ~ . + loggrow(time,
                                               model = log_growth_model))
    environment(new.cmp) <- environment()
    fit <- bru(new.cmp,
               data = data, domain = domain,
               family = family, options = list(verbose = verbose))
    print(paste("Fitted new model", n))
    n.nodes <- fit$misc$configs$nconfig
    if(!is.numeric(n.nodes)){
      print("Failed to fit, trying again")
      fit <- bru(new.cmp,
                 data = data, domain = domain,
                 family = family, options = list(verbose = verbose))
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
  log_growth_model <- define.loggrow.vary.time.model(linpoint = initial.linpoint, tmesh = tmesh, step.size = step.size, 
                                                     prior.mean = prior.mean,
                                                     prior.precision = prior.precision, priors = priors,
                                                     initial.growth = initial.growth, 
                                                     initial.carry.cap = initial.carry.cap,
                                                     initial.sigma = initial.sigma,
                                                     growth.formula = growth.formula, 
                                                     growth_cov = growth_cov,
                                                     carry.formula = carry.formula,
                                                     carry_cov = carry_cov)
  new.cmp <- update(formula, . ~ . + loggrow(time,
                                             model = log_growth_model))
  environment(new.cmp) <- environment()
  final.fit <- bru(new.cmp,
             data = data, domain = domain,
             family = family, options = list(verbose = verbose))
  return(list(fit = final.fit, n = n, linpoints = lp.mat, fit.list = fit.list))
}
