define.loggrow.joint.time.model <- function(linpoint, tmesh, step.size,
                                           prior.mean, prior.precision, nmod,
                                           priors = NULL,
                                           initial.growth = NULL, initial.carry.cap = NULL, 
                                           initial.sigma = NULL, growth.formula,
                                           carry.formula, growth_cov = NULL,
                                           carry_cov = NULL){
  #browser()
  ngrowth <- length(all.vars(growth.formula)) + attr(terms(growth.formula), "intercept") + nmod
  #assumes random effects are not in formula!
  ngrowth_cov <- length(all.vars(growth.formula))
  if(ngrowth != ngrowth_cov + nmod + 1){
    warning("Unexpected number of terms in growth formula")
  }
  #if(ngrowth_cov != ncol(growth_cov)){
  #  warning("Number of covaraites given in formula and data are not equal")
  #}
  ncarry <- length(all.vars(carry.formula)) + attr(terms(carry.formula), "intercept") + nmod
  ncarry_cov <- length(all.vars(carry.formula))
  if(ncarry != ncarry_cov + nmod + 1){
    warning("Unexpected number of terms in carry formula")
  }
  #if(ncarry_cov != ncol(carry_cov)){
  #  warning("Number of covaraites given in formula and data are not equal")
  #}
  if(is.null(initial.growth) | length(initial.growth) != ngrowth){
    initial.growth <- rep(0, ngrowth)
  }
  if(is.null(initial.carry.cap) | length(initial.carry.cap) != ncarry){
    initial.carry.cap <- rep(0, ncarry)
  }
  if(is.null(initial.sigma)){
    initial.sigma <- 0
  }
  if(length(linpoint) != nmod){
    warning("Number of linearisation points != number of models")
  }
  if(length(prior.mean) != nmod){
    warning("Number of prior means != number of models")
  }
  if(length(prior.precision) != nmod){
    warning("Number of prior precisions != number of models")
  }
  the_model <- inla.rgeneric.define(log_growth_time_joint, 
                                    linpoint = linpoint, 
                                    tmesh = tmesh, step.size = step.size, nmod = nmod,
                                    prior.mean = prior.mean, priors = priors,
                                    prior.precision = prior.precision, 
                                    initial.growth = initial.growth, 
                                    initial.carry.cap = initial.carry.cap,
                                    initial.log.sigma = initial.sigma, 
                                    ngrowth = ngrowth, ngrowth_cov = ngrowth_cov,
                                    ncarry = ncarry, ncarry_cov = ncarry_cov,
                                    growth_cov = growth_cov,
                                    carry_cov = carry_cov)
  class(the_model) <- c("time_only_joint", class(the_model))
  the_model[["tmesh"]] <- tmesh
  return(the_model)
}