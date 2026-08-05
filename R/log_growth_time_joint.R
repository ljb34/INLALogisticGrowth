log_growth_time_joint =  function(
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
    #This assumes both growth and carry have shared intercepts and individual random effects
    growth <- vector(mode = "list", length = nmod)
    carry.cap <- vector(mode = "list", length = nmod)
    for(mod in 1:nmod){
      #Shared intercepts and rand effects
      growth[[mod]] = theta[1] + theta[1+ngrowth_cov+mod]
      carry.cap[[mod]] =  theta[ngrowth+1] + theta[ngrowth + 1 + ncarry_cov + mod]
      if(ngrowth_cov > 0){
        for(i in 1:ngrowth_cov){
          growth[[mod]] <- growth[[mod]] + theta[1+i]*growth_cov[[mod]][,i]
        }
      }
      if(ncarry_cov> 0 ){
        for(i in 1:ncarry_cov){
          carry.cap[[mod]] <- carry.cap[[mod]] + theta[ngrowth + 1 + i]*carry_cov[[mod]][,i]
        }
      }
      carry.cap[[mod]]<- exp(carry.cap[[mod]])
    }
    return(list(growth = growth,
                carry.cap = carry.cap,
                sigma = exp(theta[ngrowth + ncarry + 1])))
    }
  
  graph = function() {
    return (Q())
  }
  Q = function(){
    par = interpret.theta()
    Qs <- vector(mode = "list", length = nmod)
    for(i in 1:nmod){
      Lmat = L.matrix(par$growth[[i]], par$carry.cap[[i]], step.size, linpoint[[i]], tmesh)
      #print(Lmat)
      noise.variance = Matrix::Diagonal(tmesh$n, c(prior.precision[i],rep(1/(par$sigma*step.size)**2, tmesh$n -1)))
      #print("crossprod")
      Qs[[i]] = Matrix::drop0(Matrix::crossprod(Lmat, noise.variance %*% Lmat))
    }
    output <- Matrix::bdiag(Qs)
    return(Matrix::drop0(output))
  }
  mu = function(){
    par = interpret.theta()
    output <- Matrix::Matrix(nrow = nmod*tmesh$n, ncol = 1)
    for(i in 1:nmod){
      Lmat = L.matrix(par$growth[[i]], par$carry.cap[[i]], step.size,linpoint[[i]], tmesh)
      #print(Lmat)
      r = c(prior.mean[i], r.vector(par$growth[[i]], par$carry.cap[[i]], linpoint[[i]])[-1])
      if(!is.nan(Matrix::det(Lmat))) {
        if(abs(Matrix::det(Lmat)) <= .Machine$double.eps|(is.infinite(Matrix::det(Lmat)) & !is.infinite(Matrix::det(Matrix::crossprod(Lmat,Lmat))))){ #if close to singular use
          #print(det(crossprod(Lmat,Lmat)))
          mu = Matrix::solve(Matrix::crossprod(Lmat,Lmat),Matrix::crossprod(Lmat,r)) #more stable form of solve(lmat,r)
          mu= as.vector(mu)
          #print("Trick version")
        }else{
          mu = Matrix::solve(Lmat,r)
          #print("Default Solve")
        }}else{
          #print("There's some NaNs going on?")
          mu = NA
        }
      output[((i-1)*tmesh$n + 1):(i*tmesh$n), ] <- mu
    }
    return(output)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function(){#can change params to make user specified
    #print("Calcualting logprior")
    #par = interpret.theta()
    if(is.null(priors)) warning("Parameters missing for priors")
    val = 0
    #Fixed effects
    for(i in 1:(ngrowth_cov + 1)){
      val <- val + dnorm(theta[i], priors$growth[2*i-1], priors$growth[2*i], log = T)
    }
    for(i in 1:(ncarry_cov + 1)){
      val <- val + dnorm(theta[ngrowth + i], priors$cc[2*i-1], priors$cc[2*i], log = T)
    }
    val <- val + dnorm(theta[ngrowth + ncarry + 1], mean = priors$sigma[1], sd = priors$sigma[2], log = T)
    #Random effects
    val <- val + sum(dnorm(theta[(ngrowth_cov + 2):ngrowth],0,exp(theta[ngrowth+ncarry+2]), log = T))
    val <- val + sum(dnorm(theta[(ngrowth+ncarry_cov+2):(ngrowth+ncarry)],0,exp(theta[ngrowth+ncarry+3]), log = T))
    val <- val + dnorm(theta[ngrowth+ncarry+2],mean=0,sd=1,log=TRUE)
    val <- val + dnorm(theta[ngrowth+ncarry+3],mean=0,sd=1,log=TRUE)
    return(val)
  }
  initial = function(){
    if(!exists("initial.growth", inherits = TRUE)) initial.growth = c(0.5, rep(0,ngrowth -1))
    if(!exists("initial.carry.cap", inherits = TRUE)) initial.carry.cap = c(log(1000), rep(0, ncarry-1))
    if(!exists("initial.log.sigma", inherits = TRUE)) initial.log.sigma = log(5)
    if(!exists("initial.rand.sd.growth", inherits = TRUE)) initial.rand.sd.growth = 0
    if(!exists("initial.rand.sd.carry", inherits = TRUE)) initial.rand.sd.carry = 0
    return(c(initial.growth, initial.carry.cap, initial.log.sigma, 
             initial.rand.sd.growth, initial.rand.sd.carry))
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
