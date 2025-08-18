define.cgeneric.loggrow.model <- function(linpoint, smesh, tmesh, step.size,
                                 prior.mean, prior.variance,
                                 priors = NULL, grad = NULL,
                                 initial.growth = NULL, initial.carry.cap = NULL, 
                                 initial.move.const = NULL, initial.log.sigma = NULL){
  if(is.null(grad)){
    grad <- gradient_of_linpoint(linpoint, smesh, tmesh)
  }
  mag_grad_sq <- rowSums(grad*grad)
  
  if(is.null(priors)){
    priors <- list(growth = c(0,2), cc = c(log(100), 10), #fairly flat prior
                   move = c(0,10), sigma = c(0,3))
  }
  
  fem.matrices <- fmesher::fm_fem(smesh)
  CinvG <- solve(fem.matrices$c1, fem.matrices$g1)
  
  libpath <- system.file("libs", package = "INLAloggrow")
  if (Sys.info()["sysname"] == "Windows") {
    libpath <- file.path(libpath, "INLAloggrow.dll")
  } else {
    libpath <- file.path(libpath, "INLAloggrow.so")
  }
  n <- smesh$n*tmesh*n
  
  args0 <- list(model = "inla_cgeneric_loggrow_model",
                shlib = libpath,
                n = as.integer(n),
                debug = as.integer(debug))
  
  the_model <- do.call("inla.cgeneric.define",
                       c(args0,
                         list(ns = as.integer(smesh$n),
                              nt = as.integer(tmesh$n),
                              step_size = as.double(step.size),
                              linpoint = as.double(linpoint),
                              CinvG = CinvG,
                              prior_variance = Matrix(prior.variance),
                              mag_grad_sq = as.double(mag_grad_sq),
                              prior_mean = as.double(prior.mean),
                              initial_growth = as.double(initial.growth),
                              initial_carry_cap = as.double(initial.carry.cap),
                              initial_move = as.double(initial.move.const),
                              initial_sigma = as.double(initial.log.sigma),
                              pgrowth = as.double(priors$growth),
                              pcc = as.double(priors$cc),
                              pmove = as.double(priors$move),
                              psigma = as.double(priors$sigma)
                              )))
  
  class(the_model) <- c("log_growth_model", class(the_model))
  the_model[["smesh"]] <- smesh
  the_model[["tmesh"]] <- tmesh
  return(the_model)
}