log_growth_time_vary =  function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL){ 
  
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
    a.mat <- Matrix::Diagonal(nt,a)
    
    subdiag <- Matrix::bandSparse(nt,nt,-1,list(rep(-1/step.size, nt-1)))
    
    main.diag <- Matrix::Diagonal(nt, c(0, rep(1/step.size, nt-1)))
    #print(diag(main.diag + subdiag + a.mat))
    return(Matrix::drop0(main.diag + subdiag + a.mat))
  }
  r.vector <- function(growth,carry.cap,linpoint){
    
    return(growth*(1-exp(linpoint))/carry.cap+linpoint*(1/carry.cap)*exp(linpoint))
  }
  interpret.theta = function() {
    #print(length(theta))
    growth = theta[1]
    carry.cap = theta[ngrowth+1]
    if(ngrowth > 1){
    for(i in 2:ngrowth){
      growth <- growth + theta[i]*growth_cov[,i]
    }
    }
    if(ncarry>1 ){
    for(i in 2:ncarry){
      carry.cap <- carry.cap + theta[ngrowth+i]*carry_cov[,i]
    }
    }
    return(list(growth = growth,
                carry.cap = exp(carry.cap),
                sigma = exp(theta[ngrowth + ncarry + 1])))
  }
  
  graph = function() {
    return (Q())
  }
  Q = function(){
    #print("Q being calculated")
    par = interpret.theta()
    #print(par)
    Lmat = L.matrix(par$growth, par$carry.cap, step.size, linpoint, tmesh)
    #print(Lmat)
    noise.variance = Matrix::Diagonal(tmesh$n, c(prior.precision,rep(1/(par$sigma*step.size)**2, tmesh$n -1)))
    #print("crossprod")
    output = Matrix::crossprod(Lmat, noise.variance %*% Lmat)
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
    if(!is.nan(Matrix::det(Lmat))) {
      if(abs(Matrix::det(Lmat)) <= .Machine$double.eps|(is.infinite(Matrix::det(Lmat)) & !is.infinite(Matrix::det(Matrix::crossprod(Lmat,Lmat))))){ #if close to singular use
        #print(det(crossprod(Lmat,Lmat)))
        mu = Matrix::solve(Matrix::crossprod(Lmat,Lmat),Matrix::crossprod(Lmat,r)) #more stable form of solve(lmat,r)
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
    if(is.null(priors)) warning("Parameters missing for priors")
    val = 0
    for(i in 1:ngrowth){
      val <- val + dnorm(theta[i], priors$growth[2*i-1], priors$growth[2*i], log = T)
    }
    for(i in 1:ncarry){
      val <- val + dnorm(theta[ngrowth + i], priors$cc[2*i-1], priors$cc[2*i], log = T)
    }
    val <- val + dnorm(theta[ngrowth + ncarry + 1], mean = priors$sigma[1], sd = priors$sigma[2], log = T)
    return(val)
  }
  initial = function(){
    if(!exists("initial.growth", inherits = TRUE)) initial.growth = c(0.5, rep(0,ngrowth -1))
    if(!exists("initial.carry.cap", inherits = TRUE)) initial.carry.cap = c(log(1000), rep(0, ncarry-1))
    if(!exists("initial.log.sigma", inherits = TRUE)) initial.log.sigma = log(5)
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
