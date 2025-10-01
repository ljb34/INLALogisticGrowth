#'@name define.loggrow.time.model
#'@title Define non spatial logistic growth model
#'@description
#'Defines latent model and mapper for log_growth_time
#'@param linpoint Linearisation point
#' @param tmesh mesh over equally spaced time pointes created with fm_mesh_1d
#'@param step.size difference between each of the time points
#'@param prior.mean estimated mean for first year of data. Helper function for calculating coming soon
#'@param prior.precision uncertainty for estimated mean of first year of data. Helper function for calculating coming soon
#'@param priors named list of prior parameters, named \code{cc} (carrying capacity), \code{growth}, \code{sigma}. 
#' Each is a two element vector containing the mean and variance for each parameter. 
#' @param initial.growth,initial.carry.cap,initial.log.sigma Starting values for the \emph{log} growth, 
#' \emph{log} carrying capacity and \emph{log} standard deviation
#'@returns INLA rgeneric model
#'@export


define.loggrow.time.model <- function(linpoint, tmesh, step.size,
                                              prior.mean, prior.precision, priors = NULL,
                                              initial.growth = NULL, initial.carry.cap = NULL, 
                                              initial.log.sigma = NULL){
  
  the_model <- inla.rgeneric.define(log_growth_time, 
                                    linpoint = linpoint, 
                                    tmesh = tmesh, step.size = step.size, 
                                    prior.mean = prior.mean, priors = priors,
                                    prior.precision = prior.precision, 
                                    initial.growth = initial.growth, 
                                    initial.carry.cap = initial.carry.cap,
                                    initial.log.sigma = initial.log.sigma)
  class(the_model) <- c("time_only_model", class(the_model))
  the_model[["tmesh"]] <- tmesh
  return(the_model)
}

#'@name bru_get_mapper.time_only_model
#' @title Mapper function for internal use
#' @export
bru_get_mapper.time_only_model <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  inlabru::bru_mapper(model[["tmesh"]], indexed = TRUE)
}
