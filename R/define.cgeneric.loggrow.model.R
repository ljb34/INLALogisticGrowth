#'@name define.loggrow.model
#'@title Define logistic growth model in cgeneric
#'@description
#'Defines latent model and mapper for log.growth.regeneric
#'@param linpoint Linearisation point
#'@param smesh spatial mesh created with fm_mesh_2d_inla
#' @param tmesh mesh over equally spaced time pointes created with fm_mesh_1d
#'@param step.size difference between each of the time points
#'@param prior.mean estimated mean for first year of data. Helper function for calculating coming soon
#'@param prior.precision uncertainty for estimated mean of first year of data. Helper function for calculating coming soon
#'@param priors named list of prior parameters, named \code{cc} (carrying capacity), \code{growth}, \code{move}, \code{sigma}. 
#' Each is a two element vector containing the mean and variance for each parameter. 
#' @param initial.growth,initial.carry.cap,initial.move.const,initial.log.sigma Starting values for the growth, 
#' \emph{log} carrying capacity, movement constant and \emph{log} standard deviation
#'@returns INLA cgeneric model packageCheck
#'
#'@importFrom INLAtools packageCheck
#'@export
define.cgeneric.loggrow.model <- function(linpoint, smesh, tmesh, step.size,
                                 prior.mean, prior.precision,
                                 priors = NULL, grad = NULL,
                                 initial.growth = NULL, initial.carry.cap = NULL, 
                                 initial.move.const = NULL, initial.log.sigma = NULL, debug = NULL){
  if(is.null(grad)){
    grad <- gradient_of_linpoint(linpoint, smesh, tmesh)
  }
  mag_grad_sq <- rowSums(grad*grad)
  
  if(is.null(priors)){
    priors <- list(growth = c(0,2), cc = c(log(100), 10), #fairly flat prior
                   move = c(0,10), sigma = c(0,3))
  }
  
  fem.matrices <- fmesher::fm_fem(smesh)
  CinvG <- Matrix::solve(fem.matrices$c1, fem.matrices$g1)
  #browser()
  INLAversion <- INLAtools::packageCheck(
    name = "INLA",
    minimum_version = "23.08.16",
    quietly = TRUE
  )
  
  if (INLAversion <= "25.02.10") {
    ## Old style: use INLA's external.lib helper
    libpath <- INLA::inla.external.lib("INLAloggrowth")
    hasverbose <- TRUE
  } else {
    ## Newer versions: package-installed shared library
    libpath <- system.file("libs", package = "INLAloggrowth")
    
    if (Sys.info()["sysname"] == "Windows") {
      libpath <- file.path(libpath, "INLAloggrowth.dll")
      if(!file.exists(libpath)){
        libpath <- system.file("libs", package = "INLAloggrowth")
        libpath <- file.path(libpath, "x64/INLAloggrowth.dll")
      }
    } else {
      libpath <- file.path(libpath, "INLAloggrowth.so")
    }
    hasverbose <- FALSE
  }
  stopifnot(file.exists(libpath))
  n <- smesh$n*tmesh$n
  if(is.null(debug)) debug = 0
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
                              mag_grad_sq = as.double(mag_grad_sq),
                              prior_mean = as.double(prior.mean),
                              initial_growth = as.double(initial.growth),
                              initial_carry_cap = as.double(initial.carry.cap),
                              initial_move = as.double(initial.move.const),
                              initial_sigma = as.double(initial.log.sigma),
                              pgrowth = as.double(priors$growth),
                              pcc = as.double(priors$cc),
                              pmove = as.double(priors$move),
                              psigma = as.double(priors$sigma),
                              CinvG = CinvG,
                              prior_precision = Matrix::Matrix(prior.precision, sparse = T)
                              )))
  
  class(the_model) <- c("log_growth_model", class(the_model))
  the_model[["smesh"]] <- smesh
  the_model[["tmesh"]] <- tmesh
  return(the_model)
}