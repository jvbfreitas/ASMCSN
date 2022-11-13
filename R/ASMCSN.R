print.ASMCSN = function(x, digits = NULL, quote = FALSE, prefix = "", ...){
  if(is.null(digits)){
    digits = options()$digits
  }else{
    options(digits = digits)
  }
  cat("\n", x$title)
  cat("n"," Model:\n")
  cat(" Mixed measure: ", x$mixed, "\n")
  #cat(" n:", x$n, "\n")
  cat("\nCall:\n")
  dput(x$call)                        #       cat("\nTerms:\n")
  ###        ys = matrix(rep(as.matrix(x$id, ncol = 1), 5), ncol = 5)
  #ys = matrix(rep(matrix(x$id, ncol = 1), 5), ncol = 5)
  #ys[, 2] = x$y
  #ys[, 3] = x$linear.predictors
  #ys[, 4] = x$fitted.values
  #ys[, 5] = x$residuals
  #dimnames(ys) = list(1:length(x$y), c("ID", "Y", "LP", "fitted",
  #                                      "Residual"))
  #       cat("\nFitted Values:\n")
  cat("\nNumber of observations : ", x$n, "\n")
  #cat("\nMaximum cluster size   : ", x$max.id, "\n")
  cat("\n\nParametric regression coefficients:\n")
  print(t(x$coefficients$beta), digits = digits)
  cat("\nEstimated variance parameter:\n")
  print(x$coefficients$sigma2, digits = digits)
  cat("\nEstimated skewness parameter (gamma):\n")
  print(x$coefficients$gamma, digits = digits)
  cat("\nEstimated skewness parameter (delta):\n")
  print(x$coefficients$delta, digits = digits)
  cat("\nEstimated shape parameter:\n")
  print(t(x$coefficients$nu), digits = digits)
  cat("\nNumber of Iterations: ", x$iterations)
  invisible(x)
}
summary.ASMCSN = function(object,...){
  value = list()
  class(value) = "summary.ASMCSN"
  value$call = object$call
  value$varp = object$coefficients$sigma2
  value$skew = object$coefficients$gamma
  value$shape = object$coefficients$nu
  value$alpha = object$coefficients$alpha
  value$n = object$n
  valorp = 0
  wald = 0
  pc = length(object$coefficients$beta) + 2
  cce = c(object$coefficients$beta, object$coefficients$sigma2,
          object$coefficients$delta, Reduce("c",object$coefficients$theta))
  for(i in 1:pc){
    wald[i] = wald.test(object$vcov, cce, Terms = i)$result$chi2[1]
    valorp[i] = wald.test(object$vcov, cce, Terms = i)$result$chi2[3]
  }
  value$coefficients = data.frame(cce[1:pc], object$se[1:pc], wald, valorp)
  colnames(value$coefficients) = c("Estimate","Std.err", "Wald", "Pr(>|W|)")
  for(i in 1:length(object$coefficients$beta)){
    rownames(value$coefficients)[i] = paste0("beta[",i-1,"]")
  }
  rownames(value$coefficients)[length(object$coefficients$beta)+1] = "sigma^2"
  rownames(value$coefficients)[length(object$coefficients$beta)+2] = "delta"
  return(value)
}
print.summary.ASMCSN = function(x, digits = NULL, quote = FALSE, prefix = "",
                                ...){
  if(is.null(digits)){
    digits = options()$digits
  }else{
    options(digits = digits)
  }
  cat("\nCall:\n")
  dput(x$call)                        #       cat("\nTerms:\n")
  cat("\nCoefficients:\n")
  printCoefmat(as.matrix(x$coefficients), digits = digits)
  cat("\nSample size : ", x$n, "\n")
  cat("\nEstimated shape parameter:\n")
  print(t(x$shape), digits = digits)
  cat("\nSmoothing parameters:\n")
  print(t(x$alpha), digits = digits)
  invisible(x)
}
ASMCSN = function(formula.pp, formula.npp, data, mixed = "1", k = 1/3,
                  iter.max = 300, alpha.fix = F, alpha0, nknot, a0, a1){

  namesmix = c("1", "gamma1p", "gamma2p", "beta", "binary", "betaprime",
               "birnbaum-saunders", "generalizedgamma")
  if(all(namesmix != mixed)){
    stop("The mixed measure is not defined")
  }
  if(k>1|k<0){
    stop("k must be between 0 and 1")
  }
  # if(alpha.fix==T & (alpha0 <= 0)){
  #   stop("alpha must be positive")
  # }
  formula1 = as.formula(formula.pp)
  nformula = all.vars(formula.pp)
  n2formula = all.vars(formula.npp)
  fnames = 0
  jaux = 1
  listaux = list(NULL)
  for(i in 1:length(nformula)){
    if(is.factor(data[,nformula[i]])){
      fnames[jaux] = nformula[i]
      listaux[[jaux]] = "contr.sum"
      jaux = jaux+1
    }
  }
  call = match.call()
  if(jaux>1){
    names(listaux) = fnames
    X = as.matrix(model.matrix(formula.pp, data = data, contrasts = listaux))
  }else{
    X = as.matrix(model.matrix(formula.pp, data = data))
  }
  X = as.matrix(X)
  p = ncol(X)
  y = model.frame(formula.pp, data = data)[,1]
  n = length(y)
  N = list(NULL)
  K = list(NULL)
  for(i in 1:length(n2formula)){
    t2 = ncs(data[,n2formula[i]], nknots = nknot, lambda=1)
    tr = as.vector(attr(t2, "knots"))
    K[[i]] = attr(t2, "K")
    N[[i]] = attr(t2, "N")
  }
  if(alpha.fix==F){
    alpha0 = rep(1,length(K))
    a0 = rep(0.01,length(K))
    a1 = rep(5000,length(K))
  }
  resinit = csnm(y, X, N, K, alpha0 = alpha0, nknot = nknot, tol = tol,
                 iter.max = iter.max, k = k, alpha.fix = alpha.fix, a0 = a0,
                 a1 = a1)
  mu = as.vector(resinit$fitted.values)
  sigma2 = as.vector(resinit$sigma2)
  gama = as.vector(resinit$gama)
  alpha = as.vector(resinit$alpha)
  theta = resinit$theta
  beta = as.vector(resinit$beta)
  if(mixed=="beta"){
    res = SAEM.CSS(y,mu,theta,sigma2,gama,n,p,alpha,K,N,X,k,iter.max,alpha.fix,a0,a1)
  }
  if(mixed=="gamma1p"){
    res = SAEM.CST(y,mu,theta,sigma2,gama,n,p,alpha,K,N,X,k,iter.max,alpha.fix,a0,a1)
  }
  if(mixed=="gamma2p"){
    res = SAEM.CSGT(y,mu,theta,gama,n,p,alpha,K,N,X,k,iter.max,alpha.fix,a0,a1)
  }
  if(mixed=="binary"){
    res = SAEM.CSCN(y,mu,theta,sigma2,gama,n,p,alpha,K,N,X,k,iter.max,alpha.fix,a0,a1)
  }
  if(mixed=="betaprime"){
    res = SAEM.CSBPN(y,mu,theta,sigma2,gama,n,p,alpha,K,N,X,k,iter.max,alpha.fix,a0,a1)
  }
  if(mixed=="birnbaum-saunders"){
    res = SAEM.CSBSN(y,mu,theta,sigma2,gama,n,p,alpha,K,N,X,k,iter.max,alpha.fix,a0,a1)
  }
  if(mixed=="generalizedgamma"){
    res = SAEM.CSGGN(y,mu,theta,gama,n,p,alpha,K,N,X,k,iter.max,alpha.fix,a0,a1)
  }
  if(mixed != "1"){
    mu = res$mu
    sigma2 = res$sigma2
    gama = res$gama
    alpha = res$alpha
    theta = res$theta
    beta = res$beta
    S1n = res$S1n
    S2n = res$S2n
    S3n = res$S3n
    nu = res$nu
    betao = res$betao
  }else{
    S1n = rep(1,n)
    S2n = resinit$comp$S2n
    S3n = resinit$comp$S3n
  }
  g13 = nthroot(gama,3)
  g23 = g13^2
  sigma = sqrt(sigma2)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  del = lambda/sqrt(1+lambda^2)
  Delta = sigma*del/sqrt(1-b^2*del^2)
  tau = sigma2*(1-del^2)/(1-b^2*del^2)
  A2 = -del*sqrt(1-b^2*del^2)/(2*sigma^3*(1-del^2))
  A3 = (sigma/sqrt(1-b^2*del^2))*(1+del^2*b^2/(1-b^2*del^2))
  A4 = (2*del*sigma2/(1-b^2*del^2))*((1-del^2)*b^2/(1-b^2*del^2)-1)
  A5 = (tau*A3-Delta*A4)/(tau^2)
  A6 = (2*tau*Delta*A3-Delta^2*A4)/(tau^2)
  if(mixed == "gamma1p"){
    vary = sigma2*nu/(nu-2)
  }
  if(mixed == "gamma2p"){
    vary = sigma2*nu[2]/(nu[1]-2)
  }
  if(mixed == "binary"){
    vary = sigma2*(nu[1]+nu[2]*(1-nu[1]))/nu[2]
  }
  if(mixed == "beta"){
    vary = sigma2*nu/(nu-1)
  }
  if(mixed == "betaprime"){
    vary = sigma2*nu[2]/(nu[1]-1)
  }
  if(mixed == "birnbaum-saunders"){
    vary = sigma2*(nu[1]^2+2)/(2*nu[2])
  }
  if(mixed == "generalizedgamma"){
    vary = sigma2*gamma(nu[1]-1/nu[2])/(nu[3]*gamma(nu[1]))
  }
  if(mixed == "1"){
    vary = sigma2
  }
  sXi = list(NULL)
  stheta = list(NULL)
  for(i in 1:n){
    sbeta = -(1/tau)*(t(X[i,])*Delta*S3n[i]-t(X[i,])*(y[i]-mu[i])*S1n[i])
    for(j in 1:length(N)){
      stheta[[j]] = -(1/tau)*((N[[j]][i,])*Delta*S3n[i]-t(N[[j]][i,])*(y[i]-mu[i])*S1n[i]) - t(.5*(alpha[[j]]*K[[j]]%*%(theta[[j]])))
    }
    ss2 = -(1/(2*sigma2))+.5*(1/(tau*sigma2))*(y[i]-mu[i])^2*S1n[i]+ 2*(y[i]-mu[i])*S3n[i]*A2
    sdelta = -1/(2*tau)*A4 -(.5)*(-(1/(tau^2))*A4*(y[i]-mu[i])^2*S1n[i]-2*A5*(y[i]-mu[i])*S3n[i]+A6*S2n[i])
    sXi[[i]] = c(sbeta,ss2,sdelta,Reduce("c", stheta))%*%t(c(sbeta,ss2,sdelta,Reduce("c", stheta)))
  }
  # sXi = list(NULL)
  # for(i in 1:n){
  #   sbeta = -(1/tau)*(t(x[i,])*Delta*S3n[i]-t(x[i,])*(y[i]-mu[i])*S1n[i])
  #   stheta = -(1/tau)*((N[i,])*Delta*S3n[i]-t(N[i,])*(y[i]-mu[i])*S1n[i]) - t(.5*(alpha*K%*%(theta)))
  #   sdelta = -1/(2*tau)*A4 -(.5)*(-(1/(tau^2))*A4*(y[i]-mu[i])^2*S1n[i]-2*A5*(y[i]-mu[i])*S3n[i]+A6*S2n[i])
  #   sXi[[i]] = c(sbeta,stheta,ss2,sdelta)%*%t(c(sbeta,stheta,ss2,sdelta))
  # }
  Ie = Reduce("+", sXi)
  avcov = ginv(Ie)
  vare = sqrt(diag(avcov))
  fit = list()
  attr(fit, "class") = c("ASMCSN")
  if(mixed == "gamma1p"){
    fit$title = "CST:  CENTERED SKEW-T MODEL"
    fit$type = "CST"
  }
  if(mixed == "gamma2p"){
    fit$title = "CSGT:  CENTERED SKEW GENERALIZED T MODEL"
    fit$type = "CSGT"
  }
  if(mixed == "binary"){
    fit$title = "CSCN:  CENTERED SKEW CONTAMINATED NORMAL MODEL"
    fit$type = "CSCN"
  }
  if(mixed == "beta"){
    fit$title = "CSS:  CENTERED SKEW SLASH MODEL"
    fit$type = "CSS"
  }
  if(mixed == "betaprime"){
    fit$title = "CSBPN:  CENTERED SKEW BETA PRIME NORMAL MODEL"
    fit$type = "CSBPN"
  }
  if(mixed == "birnbaum-saunders"){
    fit$title = "CSBSN:  CENTERED SKEW BIRNBAUM-SAUNDERS NORMAL MODEL"
    fit$type = "CSBSN"
  }
  if(mixed == "generalizedgamma"){
    fit$title = "CSGGN:  CENTERED SKEW GENERALIZED GAMMA NORMAL MODEL"
    fit$type = "CSGGN"
  }
  if(mixed == "1"){
    fit$title = "CSN:  CENTERED SKEW-NORMAL MODEL"
    fit$type = "CSN"
  }
  fit$call = call
  fit$formula.pp = formula.pp
  fit$formula.npp = formula.npp
  fit$n = n
  fit$mixed = mixed
  fit$iterations = iter.max
  fit$coefficients$beta = beta
  fit$coefficients$theta = theta
  fit$coefficients$gamma = gama
  fit$coefficients$delta = del
  fit$coefficients$sigma2 = sigma2
  fit$coefficients$alpha = alpha
  fit$coefficients$tau = tau
  if(mixed != 1){
    fit$coefficients$nu = nu
  }
  fit$fitted.values = mu
  r = (y-mu)/sigma
  fit$residuals = r
  fit$mahalanobis = (y-mu)^2/vary
  fit$y = as.vector(y)
  fit$vcov = avcov
  fit$se = vare
  fit$comp$X = X
  fit$comp$N = N
  fit$comp$K = K
  fit$expectations$S1n = S1n
  fit$expectations$S2n = S2n
  fit$expectations$S3n = S3n
  if(mixed != "1"){
    fit$trajectory$betao = res$betao
    if(mixed == "gamma2p"|mixed == "generalizedgamma"){
      fit$trajectory$sigma2o = rep(1,iter.max)
    }else{
      fit$trajectory$sigma2o = res$sigma2o
    }
    fit$trajectory$gammao = res$gamao
    fit$trajectory$nuo = res$nuo
  }else{
    fit$trajectory$betao = resinit$betao
    fit$trajectory$sigma2o = resinit$sigma2o
    fit$trajectory$gammao = resinit$gamao
  }
  if(mixed == "gamma1p"){
    lp = sum(log(dCST(y, mu, sigma2, gama, nu)))
  }
  if(mixed == "gamma2p"){
    lp = sum(log(dCSGT(y, mu, gama, nu)))
  }
  if(mixed == "binary"){
    lp = sum(log(dCSCN(y, mu, sigma2, gama, nu)))
  }
  if(mixed == "beta"){
    lp = sum(log(dCSS(y, mu, sigma2, gama, nu)))
  }
  if(mixed == "betaprime"){
    lp = sum(log(dCSBPN(y, mu, sigma2, gama, nu)))
  }
  if(mixed == "birnbaum-saunders"){
    lp = sum(log(dCSBSN(y, mu, sigma2, gama, nu)))
  }
  if(mixed == "generalizedgamma"){
    lp = sum(log(dCSGGN(y, mu, gama, nu)))
  }
  if(mixed == "1"){
    lp = sum(log(dCSN(y, mu, sigma2, gama)))
  }
  df = rep(0,length(K))
  fit$logver = lp
  for(i in 1:length(alpha)){
    q=ncol(K[[i]])
    L=tau*alpha[i]*K[[i]]
    aux = N[[i]]%*%solve(t(N[[i]])%*%diag(S1n)%*%N[[i]]+L)%*%t(N[[i]])
    df[i]=sum(diag(aux))
    lp = lp - 0.5*alpha[i]*(t(theta[[i]])%*%K[[i]]%*%theta[[i]])
  }
  fit$df = df
  return(fit)
}

