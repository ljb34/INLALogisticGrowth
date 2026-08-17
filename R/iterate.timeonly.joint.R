iterate.timeonly.joint<- function(formula = ~-1, data, tmesh, step.size, prior.mean, 
                                  nmod, group_identifier,
                                  prior.precision, max.iter = 100,gamma = 0.75,stop.crit = 0.05,
                                  priors = NULL, growth.formula = ~1, 
                                  carry.formula = ~1, covariates = NULL,
                                  initial.linpoint = NULL, initial.growth=1, 
                                  initial.carry.cap=0.05, initial.sigma = log(1.5), 
                                  options = list(verbose = F), family = "gaussian", domain = NULL){

  if(is.null(initial.linpoint)){
    initial.linpoint <- list()
  }
  if(length(initial.linpoint) != nmod){
    for(i in 1:nmod){
      initial.linpoint[[i]] <- log(logit.nest(exp(prior.mean[i]), initial.growth[1], exp(initial.carry.cap[1]), tmesh$n)$x)
    }
  }
  fit.list <- list()
  if(is.null(domain)){
    domain = list(time = seq(min(data$time), max(data$time), by = step.size))
  }
  
  if(!(group_identifier %in% names(data))){
    warning("Missing group identifier in data")
  }
  #Arrange data
  data_arrange <- dplyr::arrange(data, group_identifier)
  unique_ids <- unique(data[[group_identifier]])
  if(nmod != length(unique_ids)){
    warning("Number of models and number of groups not equal")
  }
  #Covariates
  vars_growth <- all.vars(growth.formula)
  vars_carry  <- all.vars(carry.formula)
  #Assume there is always an intercept in this model, so no need for a column of 1s
  if(length(vars_growth) >= 1){
    growth_cov <- vector(mode = "list", nmod)
    for(mod in 1:nmod){
      growth_cov[[mod]] <- matrix(nrow = tmesh$n, ncol = length(vars_growth))
      for(i in 1:length(vars_growth)){
        if(vars_growth[i] %in% names(data)){
          nvals <- sum(data[[group_identifier]] == unique_ids[mod])
          if(nvals == tmesh$n){
            growth_cov[[mod]][,i] <- data[data[[group_identifier]] == unique_ids[mod],][[vars_growth[i]]]
          } else{
            growth_cov[[mod]][,i] <- c(data[data[[group_identifier]] == unique_ids[mod],][[vars_growth[i]]], 
                                      rep(data[data[[group_identifier]] == unique_ids[mod],][[vars_growth[i]]][nvals], tmesh$n - nvals))
          }
        } else {
          warning(paste("Couldn't find covariate in dataframe", vars_growth[i]))
        }
      }
    }
    }else{
    growth_cov = NULL
  }
  if(length(vars_carry) >= 1){
    carry_cov <- vector(mode = "list", nmod)
    for(mod in 1:nmod){
      carry_cov[[mod]] <- matrix(nrow = tmesh$n, ncol = length(vars_carry))
      for(i in 1:length(vars_carry)){
        if(vars_carry[i] %in% names(data)){
          nvals <- sum(data[[group_identifier]] == unique_ids[mod])
          if(nvals == tmesh$n){
            carry_cov[[mod]][,i] <- data[data[[group_identifier]] == unique_ids[mod],][[vars_carry[i]]]
          } else{
            carry_cov[[mod]][,i] <- c(data[data[[group_identifier]] == unique_ids[mod],][[vars_carry[i]]], 
                                      rep(data[data[[group_identifier]] == unique_ids[mod],][[vars_carry[i]]][nvals], tmesh$n - nvals))
          }
        } else {
          warning(paste("Couldn't find covariate in dataframe", vars_carry[i]))
        }
      }
    }
  }else{
    carry_cov = NULL
  }
  #Set up initial model
  log_growth_model <- define.loggrow.joint.time.model(linpoint = initial.linpoint, tmesh = tmesh, step.size = step.size, 
                                                     prior.mean = prior.mean,nmod = nmod,
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
  #browser()
  fit <- bru(new.cmp,
             data = data_arrange, domain = domain,
             family = family, options = options)
  fit.list[[1]]<-fit
  print("First fitting finished")
  n.nodes <- fit$misc$configs$nconfig
  nodes <- data.frame(log.prob=rep(NA,n.nodes))
  mat_list <- list()
  mean_list <- list()
  for(i in 1:n.nodes){
    nodes[i,]<- c(fit$misc$configs$config[[i]]$log.posterior)
    mat_list[[i]] <- fit$misc$configs$config[[i]]$Q[1:(tmesh$n*nmod), 1:(tmesh$n*nmod)]
    mean_list[[i]] <- fit$misc$configs$config[[i]]$improved.mean[1:(tmesh$n*nmod)]
  }
  nodes <- mutate(nodes, weight = exp(log.prob)) %>%
    mutate(weight.prob = weight/sum(weight))
  #Old rule- in theory faster but gives some extreme changes
  #P <- Reduce("+", Map(function(m, w) m * w, mat_list, nodes$weight.prob))
  #weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
  #b <- Reduce("+", Map(function(m,w) m%*%w, mat_list,weighted.means))
  #new.linpoint <- (1-gamma)*initial.linpoint +gamma*solve(P,b)
  #browser()
  #New update rule
  weighted.means <- Map(function(v,p) v*p, mean_list, nodes$weight.prob)
  new.mean <- Reduce("+", weighted.means)
  #print(new.mean)
  new.linpoint <- list()
  for(i in 1:nmod){
    new.linpoint[[i]] <- (1-gamma)*initial.linpoint[[i]] +gamma*new.mean[(i-1)*tmesh$n+1:tmesh$n]
  }
  #Check that this linpoint isn't so extreme that it will cause issues
  #plot(new.linpoint)
  lp.mat <- cbind(unlist(initial.linpoint),unlist(new.linpoint))
  n <- 2
  #print(fit$summary.hyperpar$mean)
  
  #Iterate the updates
  while(n < max.iter & mean(abs(lp.mat[,n]-lp.mat[,n-1]))>stop.crit){
    log_growth_model <- define.loggrow.joint.time.model(linpoint = new.linpoint, tmesh = tmesh, step.size = step.size, 
                                                        prior.mean = prior.mean,nmod = nmod,
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
               data = data_arrange, domain = domain,
               family = family, options = options)
    print(paste("Fitted new model", n))
    n.nodes <- fit$misc$configs$nconfig
    if(!is.numeric(n.nodes)){
      print("Failed to fit, trying again")
      fit <- bru(new.cmp,
                 data = data_arrange, domain = domain,
                 family = family, options = options)
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
      mat_list[[i]] <- fit$misc$configs$config[[i]]$Q[1:(tmesh$n*nmod),1:(tmesh$n*nmod)]
      mean_list[[i]] <- fit$misc$configs$config[[i]]$improved.mean[1:(tmesh$n*nmod)]
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
    for(i in 1:nmod){
      new.linpoint[[i]] <- (1-gamma)*new.linpoint[[i]] +gamma*new.mean[(i-1)*tmesh$n+1:tmesh$n]
    }
    #plot(new.linpoint, main = paste("Linearisation point", n))
    lp.mat <- cbind(lp.mat,unlist(new.linpoint))
    print("Updated linpoint")
    n <- n+1
  }
  log_growth_model <- define.loggrow.joint.time.model(linpoint = new.linpoint, tmesh = tmesh, step.size = step.size, 
                                                      prior.mean = prior.mean,nmod = nmod,
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
                   data = data_arrange, domain = domain,
                   family = family, options = options)
  return(list(fit = final.fit, n = n, linpoints = lp.mat, fit.list = fit.list))
}
