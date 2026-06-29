#'Extract prior mean and variance
#'@description
#'A helper function for getting prior mean and variance from an \strong{intercept and smoothing term only} inlabru model
#'fitted to first year of data. For more complicated models (e.g. including covariates or thinning),
#'the user needs to calculate these 
#'@param fit A fitted inlabru model of class "bru". 
#'@param smesh Spatial mesh used to fit model
#'@importFrom stringr str_sub
#' @importFrom Matrix Diagonal
#'@returns named list of prior mean and variance
#'@export

get.initial.mean.var <- function(fit, smesh){
  
  index <- min(which(stringr::str_sub(rownames(fit$summary.fitted.values),8,8)!= "A"))
  initial.mean <- fit$summary.fixed$mean +fit$summary.fitted.values$mean[index-1 +1:smesh$n]
  initial.variance <- Matrix::Diagonal(mesh_obs$n, (fit$summary.fixed$sd**2)+(fit$summary.fitted.values$sd[index-1 +1:smesh$n]**2))
  return(list(prior.mean = initial.mean, prior.variance = initial.variance))
}