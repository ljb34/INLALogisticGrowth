#'rgeneric Spatial Logistic Growth
#'
#'@description
#'The rgeneric function implementing the 2D spatial logistic growth function. 
#'For use within INLA only. See inla.doc('rgeneric') for more
#'
#'@export

log_growth_rgeneric =  function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL){ 
  envir = parent.env(environment()) #gets extra parameters (linpoint etc.) from definition data

  library(Matrix)
  library(fmesher)
  a.func <- function(growth,carry.cap, linpoint){
    #print("Calcualting a")
    return(growth*exp(linpoint)/carry.cap)
  }
  #growth, carry.cap, move.const = theta params to be est
  #step.size = difference in time between lin points, known
  #linpoint = list of linearisation point vectors 
  #smesh = space mesh built with fmesher, tmesh = time mesh
  L.matrix <- function(growth,carry.cap,move.const,step.size, linpoint, smesh, tmesh){
    #print("Calcualting Lmat")
    ns <- smesh$n
    nt <- tmesh$n
    a<- a.func(growth,carry.cap, linpoint)
    a[1:ns] <- 1
    a.mat <- Matrix::Diagonal(ns*nt,a)
    subdiag <- Matrix::bandSparse(nt*ns, k = -ns, diagonals = list(rep(-1/step.size, (nt - 1)*ns)))
    fem.matrices <- fmesher::fm_fem(smesh)
    CinvG <- solve(fem.matrices$c1, fem.matrices$g1)
    main.diag <- Matrix::kronecker(Matrix::Diagonal(nt, c(0,rep(1, nt-1))), 
                                   Matrix::Diagonal(ns, 1/(step.size))+ move.const*CinvG)
    #print(diag(main.diag + subdiag + a.mat))
    return(drop0(main.diag + subdiag + a.mat, 1e-100))
  }
  r.vector <- function(growth,carry.cap,move.const,linpoint,grad){
    mag.grad.sq <- rowSums(grad*grad) #magnitude squared
    return(growth*exp(linpoint)*(linpoint-1)/carry.cap+ growth - move.const*mag.grad.sq)
  }
  interpret.theta = function() {
    return(list(growth = exp(theta[1L]),
                carry.cap = exp(theta[2L]),
                move.const = theta[3L], 
                sigma = exp(theta[4L])))
  }
  
  graph = function() {
    return (Q())
  }
  Q = function(){
    #print("Calcualting Q")
    par = interpret.theta()
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, par$move.const,step.size, linpoint, smesh, tmesh)
    noiseonly = Matrix::Diagonal(smesh$n*(tmesh$n-1), (par$sigma*step.size)**-2)
    noise.precision = Matrix::bdiag(list(prior.precision, noiseonly))
    output = Matrix::crossprod(Lmat, noise.precision %*% Lmat)
    #print(output[smesh$n:(smesh$n +10),smesh$n:(smesh$n +10)])
    return(drop0(output,1e-100))
  }
  mu = function(){
    #browser()
    #print("Calcualting mu")
    #if(class(theta)!="numeric"){
    #  theta <- initial()
    #}
    par = interpret.theta()
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, par$move.const, step.size, linpoint, smesh, tmesh)
    r = c(prior.mean, r.vector(par$growth, par$carry.cap, par$move.const, linpoint, grad)[-(1:smesh$n)])
    #print(det(Lmat))
    if(!is.nan(det(Lmat))) {
      if(abs(det(Lmat)) <= .Machine$double.eps|(is.infinite(det(Lmat)) & !is.infinite(det(crossprod(Lmat,Lmat))))){ #if close to singular use
        #print(det(crossprod(Lmat,Lmat)))
        mu = Matrix::solve(Matrix::crossprod(Lmat,Lmat),Matrix::crossprod(Lmat,r)) #more stable form of solve(lmat,r)
        mu= as.vector(mu)
        #print("Trick version")
      }else{
        mu = solve(Lmat,r)
        #print("Default Solve")
      }}else{
        print("There's some NaNs going on?")
        mu = NA
      }
    #print(mean(mu))
    return(mu)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function(){#can change params to make user specified
    #print("Calcualting logprior")
    par = interpret.theta()
    #print(par)
    if(!is.null(priors)) warning("Parameters missing for priors")
    val = dnorm(par$carry.cap, mean = priors$cc[1], sd = priors$cc[2], log = T)+
      dnorm(par$growth, mean = priors$growth[1], sd = priors$growth[2], log = T)+
      dnorm(par$move.const,mean = priors$move[1], sd = priors$move[2], log = T)+ 
      dnorm(par$sigma, mean = priors$sigma[1], sd = priors$sigma[2], log = T)
    return(val)
  }
  initial = function(){#can change params to make user specified
    if(is.null(initial.growth)) initial.growth = log(0.5)
    if(is.null(initial.carry.cap)) initial.carry.cap = log(1000)
    if(is.null(initial.move.const)) initial.move.const = 1
    if(is.null(initial.log.sigma)) initial.log.sigma = log(0.5)
    return(c(initial.growth, initial.carry.cap, 
             initial.move.const, initial.log.sigma))
  }
  quit = function() {
    return(invisible())
  }
  if (is.null(theta)) theta = initial()
  if (length(theta) == 0) theta = initial()
  val = do.call(match.arg(cmd), args = list())
  return(val)
}

