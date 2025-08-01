#' Non spatial 1 year logistic growth function
#'  Deterministic function for calculating population size after 1 time step of logistic growth
#'  @param x numeric or vector starting population
#'  @param r numeric or vector growth rate
#'  @param k numeric or vector carrying capacity
#'  @returns population after 1 time step of logistic growth
#'  @export

logit.growth <- function(x,r,k){
  x[x<0]<-0
  xnew <- x*exp(r*(1-(x/k)))
  return(xnew)
}
#'Nested logistic growth
#'Deterministic function for calculating population size over multiple time steps of growth
#'  @param x0 numeric or vector starting population
#'  @param r numeric or vector growth rate
#'  @param k numeric or vector carrying capacity
#'  @param n number of time steps to compute
#'  @returns dataframe of population size and time
#'  @export
logit.nest <- function(x0,r,k,n){
  df <- data.frame(x = x0, time = rep(1, length(x0)))
  for(i in 2:n){
    df <- rbind(df,
                data.frame(x = logit.growth(df$x[df$time == i-1], r,k), 
                           time = rep(i, length(x0))))
  }
  return(df)
}
