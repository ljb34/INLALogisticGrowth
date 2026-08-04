#'@name define.loggrow.time.model
#'@title Define non spatial logistic growth model with varying parameters
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


define.loggrow.vary.time.model <- function(linpoint, tmesh, step.size,
                                      prior.mean, prior.precision, priors = NULL,
                                      initial.growth = NULL, initial.carry.cap = NULL, 
                                      initial.sigma = NULL, growth.formula,
                                      carry.formula, growth_cov = NULL,
                                      carry_cov = NULL){
  #browser()
  ngrowth <- length(all.vars(growth.formula)) + attr(terms(growth.formula), "intercept")
  ncarry <- length(all.vars(carry.formula))+ attr(terms(carry.formula), "intercept")
  if(is.null(initial.growth) | length(initial.growth) != ngrowth){
    initial.growth <- rep(0, ngrowth)
  }
  if(is.null(initial.carry.cap) | length(initial.carry.cap) != ncarry){
    initial.carry.cap <- rep(0, ncarry)
  }
  if(is.null(initial.sigma)){
    initial.sigma <- 0
  }
  the_model <- inla.rgeneric.define(log_growth_time_vary, 
                                    linpoint = linpoint, 
                                    tmesh = tmesh, step.size = step.size, 
                                    prior.mean = prior.mean, priors = priors,
                                    prior.precision = prior.precision, 
                                    initial.growth = initial.growth, 
                                    initial.carry.cap = initial.carry.cap,
                                    initial.log.sigma = initial.sigma, 
                                    ngrowth = ngrowth,
                                    ncarry = ncarry, growth_cov = growth_cov,
                                    carry_cov = carry_cov)
  class(the_model) <- c("time_only_vary", class(the_model))
  the_model[["tmesh"]] <- tmesh
  return(the_model)
}

#'@name bru_get_mapper.time_only_model
#' @title Mapper function for internal use
#' @export
bru_get_mapper.time_only_vary <- function(model, ...) {
  stopifnot(requireNamespace("inlabru"))
  inlabru::bru_mapper(model[["tmesh"]], indexed = TRUE)
}
