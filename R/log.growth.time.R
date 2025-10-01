#'rgeneric Time only Logistic Growth (not spatial)
#'
#'@description
#'The rgeneric function implementing the time only (non-spatial) logistic growth function. 
#'For use within INLA only. See inla.doc('rgeneric') for more
#'
#'@export
#'
log_growth_time =  function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL){ 
  envir = parent.env(environment()) #gets extra parameters (linpoint etc.) from definition data
  library(Matrix)#shouldn't have to put these inside this function
  library(fmesher)#but I can't get it to work otherwise...
  
  #growth, inv.carry.cap, move.const = theta params to be est
  #step.size = difference in time between lin points, known
  #linpoint = list of linearisation point vectors, tmesh = time mesh
  
  a.func <- function(growth,carry.cap, linpoint){
    #print("a func")
    return(growth*exp(linpoint)*(1/carry.cap))
  }
  
  L.matrix <- function(growth,carry.cap,step.size, linpoint, tmesh){
    #print("Lmat start")
    nt <- tmesh$n
    a<- a.func(growth,carry.cap,linpoint)
    a[1] <- 1
    a.mat <- Diagonal(nt,a)
    
    subdiag <- bandSparse(nt,nt,-1,list(rep(-1/step.size, nt-1)))
    
    main.diag <- Diagonal(nt, c(0, rep(1/step.size, nt-1)))
    #print(diag(main.diag + subdiag + a.mat))
    return(Matrix::drop0(main.diag + subdiag + a.mat))
  }
  r.vector <- function(growth,carry.cap,linpoint){
    
    return(growth*(1-exp(linpoint))/carry.cap+linpoint*(1/carry.cap)*exp(linpoint))
  }
  interpret.theta = function() {
    #print(theta)
    return(list(growth = exp(theta[1L]),
                carry.cap = exp(theta[2L]),
                sigma = exp(theta[3L])))
  }
  
  graph = function() {
    #print("It's graph time")
    return (Q())
  }
  Q = function(){
    #print("Q being calculated")
    par = interpret.theta()
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, step.size, linpoint, tmesh)
    #print(Lmat)
    noise.variance = Diagonal(tmesh$n, c(prior.precision,rep(1/(par$sigma*step.size)**2, tmesh$n -1)))
    #print("crossprod")
    output = crossprod(Lmat, noise.variance %*% Lmat)
    #print("finished Q")
    #print(output)
    return(Matrix::drop0(output))
  }
  mu = function(){
    #print("mu being calculated")
    #if(class(theta)!="numeric"){
    #  theta <- initial()
    #}
    par = interpret.theta()
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, step.size,linpoint, tmesh)
    #print(Lmat)
    r = c(prior.mean, r.vector(par$growth, par$carry.cap, linpoint)[-1])
    #print(r)
    #print(det(Lmat))
    if(!is.nan(det(Lmat))) {
      if(abs(det(Lmat)) <= .Machine$double.eps|(is.infinite(det(Lmat)) & !is.infinite(det(crossprod(Lmat,Lmat))))){ #if close to singular use
        #print(det(crossprod(Lmat,Lmat)))
        mu = solve(crossprod(Lmat,Lmat),crossprod(Lmat,r)) #more stable form of solve(lmat,r)
        mu= as.vector(mu)
        #print("Trick version")
      }else{
        mu = solve(Lmat,r)
        #print("Default Solve")
      }}else{
        #print("There's some NaNs going on?")
        mu = NA
      }
    #print(mu)
    return(mu)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function(){#can change params to make user specified
    #print("Calcualting logprior")
    par = interpret.theta()
    if(!is.null(priors)) warning("Parameters missing for priors")
    val = dnorm(theta[2L], mean = priors$cc[1], sd = priors$cc[2], log = T)+ 
      dnorm(theta[1L], mean = priors$growth[1], sd = priors$growth[2], log = T)+
      dnorm(theta[3L], mean = priors$sigma[1], sd = priors$sigma[2], log = T)
    return(val)
  }
  initial = function(){
    if(is.null(initial.growth)) initial.growth = 0.5
    if(is.null(initial.carry.cap)) initial.carry.cap = 1000
    if(is.null(initial.log.sigma)) initial.log.sigma = log(5)
    return(c(initial.growth, initial.carry.cap, initial.log.sigma))
  }
  quit = function() {
    return(invisible())
  }
  if (is.null(theta)) theta = initial()
  if (length(theta) == 0) theta = initial()
  #if (NaN %in% theta) print(theta)
  val = do.call(match.arg(cmd), args = list())
  return(val)
}