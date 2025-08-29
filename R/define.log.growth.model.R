#'@name define.loggrow.model
#'@title Define logistic growth model
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
#'@returns INLA rgeneric model
#'@export
define.loggrow.model <- function(linpoint, smesh, tmesh, step.size,
                                 prior.mean, prior.precision,
                                priors = NULL, grad = NULL,
                                 initial.growth = NULL, initial.carry.cap = NULL, 
                                 initial.move.const = NULL, initial.log.sigma = NULL){
  if(is.null(grad)){
    grad <- gradient_of_linpoint(linpoint, smesh, tmesh)
  }
  the_model <- inla.rgeneric.define(log_growth_rgeneric, 
                                    linpoint = linpoint, 
                                    smesh = smesh,tmesh = tmesh, step.size = step.size, 
                                    prior.mean = prior.mean, priors = priors,
                                    prior.precision = prior.precision, 
                                    grad = grad,
                                    initial.growth = initial.growth, 
                                    initial.carry.cap = initial.carry.cap, 
                                    initial.move.const = initial.move.const, 
                                    initial.log.sigma = initial.log.sigma)
  class(the_model) <- c("log_growth_model", class(the_model))
  the_model[["smesh"]] <- smesh
  the_model[["tmesh"]] <- tmesh
  return(the_model)
}

#'@name bru_get_mapper.log_growth_model
#' @title Mapper function for internal use
#' @export
bru_get_mapper.log_growth_model <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  inlabru::bru_mapper_multi(list(
    space = inlabru::bru_mapper(model[["smesh"]]),
    time = inlabru::bru_mapper(model[["tmesh"]], indexed = TRUE)
  ))
}
