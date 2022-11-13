qq.ASMCSN = function(obj, confi = 0.95, repl = 100){
  #opcaoEnvelope = "Normal"
  res = obj$residuals
  n = length(res)
  rres = matrix(0, nrow = n, ncol = repl)
  sres = sort(res)
  mixed = obj$mixed
  k = 0
  if(mixed == "gamma1p"){
    gama = obj$coefficients$gamma
    nu = obj$coefficients$nu
    while(k<repl){
      k = k+1
      aux20=rCST(n,0,1,gama,nu)
      rres[,k] = aux20
    }
  }
  if(mixed == "gamma2p"){
    gama = obj$coefficients$gamma
    nu = obj$coefficients$nu
    while(k<repl){
      k = k+1
      aux20=rCSGT(n,0,gama,nu)
      rres[,k] = aux20
    }
  }
  if(mixed == "binary"){
    gama = obj$coefficients$gamma
    nu = obj$coefficients$nu
    while(k<repl){
      k = k+1
      aux20=rCSCN(n,0,1,gama,nu)
      rres[,k] = aux20
    }
  }
  if(mixed == "beta"){
    gama = obj$coefficients$gamma
    nu = obj$coefficients$nu
    while(k<repl){
      k = k+1
      aux20=rCSS(n,0,1,gama,nu)
      rres[,k] = aux20
    }
  }
  if(mixed == "betaprime"){
    gama = obj$coefficients$gamma
    nu = obj$coefficients$nu
    while(k<repl){
      k = k+1
      aux20=rCSBPN(n,0,1,gama,nu)
      rres[,k] = aux20
    }
  }
  if(mixed == "birnbaum-saunders"){
    gama = obj$coefficients$gamma
    nu = obj$coefficients$nu
    while(k<repl){
      k = k+1
      aux20=rCSBSN(n,0,1,gama,nu)
      rres[,k] = aux20
    }
  }
  if(mixed == "generalizedgamma"){
    gama = obj$coefficients$gamma
    nu = obj$coefficients$nu
    while(k<repl){
      k = k+1
      aux20=rCSGGN(n,0,gama,nu)
      rres[,k] = aux20
    }
  }
  if(mixed == "1"){
    gama = obj$coefficients$gamma
    while(k<repl){
      k = k+1
      aux20=rCSN(n,0,1,gama)
      rres[,k] = aux20
    }
  }

  srres = matrix(0, nrow = n, ncol = repl)

  for(k in 1:repl) {srres[,k]=sort(rres[,k])}
  descq = matrix(0, nrow = n, ncol = 3)
  for(k in 1:n){
    descq[k,1] = quantile(srres[k,],probs = 1-confi)
    descq[k,2] = median(srres[k,])
    descq[k,3] = quantile(srres[k,],probs = confi)
  }
  if(mixed == "gamma1p"){
    aux50 = qCST(ppoints(n),0,1,gama,nu)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Standard Residuals", x = "CST quantiles")+
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(mixed == "gamma2p"){
    aux50 = qCSGT(ppoints(n),0,gama,nu)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Standard Residuals", x = "CSGT quantiles")+
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(mixed == "binary"){
    aux50 = qCSCN(ppoints(n),0,1,gama,nu)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Standard Residuals", x = "CSCN quantiles")+
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(mixed == "beta"){
    aux50 = qCSS(ppoints(n),0,1,gama,nu)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Standard Residuals", x = "CSS quantiles")+
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(mixed == "betaprime"){
    aux50 = qCSBPN(ppoints(n),0,1,gama,nu)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Standard Residuals", x = "CSBPN quantiles")+
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(mixed == "birnbaum-saunders"){
    aux50 = qCSBSN(ppoints(n),0,1,gama,nu)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Standard Residuals", x = "CSBSN quantiles")+
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(mixed == "generalizedgamma"){
    aux50 = qCSGGN(ppoints(n),0,gama,nu)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Standard Residuals", x = "CSGGN quantiles")+
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  if(mixed == "1"){
    aux50 = qCSN(ppoints(n),0,1,gama)
    ptg = ggplot() + geom_point(aes(x = aux50, y = sres)) +
      geom_line(aes(x = aux50, y = descq[,1])) +
      geom_line(aes(x = aux50, y = descq[,3])) +
      geom_line(aes(x = aux50, y = descq[,2]), linetype = "dashed") +
      labs(y = "Standard Residuals", x = "CSN quantiles")+
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text=element_text(size=15)
      )
  }
  print(ptg)
}

diag.ASMCSN = function(obj, type = c("Gcooks", "Ghat")){
  X = obj$comp$X
  N = obj$comp$N
  y = obj$y
  K = obj$comp$K
  n = nrow(X)

  alpha = obj$coefficients$alpha
  mu = obj$fitted.values
  beta = obj$coefficients$beta
  theta = obj$coefficients$theta
  sigma2 = obj$coefficients$sigma2
  gama = obj$coefficients$gamma
  nu = obj$coefficients$nu
  g13 = nthroot(gama,3)
  g23 = g13^2
  sigma = sqrt(sigma2)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  del = lambda/sqrt(1+lambda^2)
  Delta = sigma*del/sqrt(1-b^2*del^2)
  tau = sigma2*(1-del^2)/(1-b^2*del^2)
  S1n = obj$expectations$S1n
  S2n = obj$expectations$S2n
  S3n = obj$expectations$S3n
  q11 = -(1/tau)*t(X)%*%diag(S1n)%*%X
  q12 = list(NULL)
  q21 = list(NULL)
  for(i in 1:length(N)){
    q12[[i]] = -(1/tau)*t(X)%*%diag(S1n)%*%N[[i]]
    q21[[i]] = t(q12[[i]])
  }
  q13 = -(1/(tau^2))*(t(X)%*%diag(S1n)%*%(y-mu)-Delta*t(X)%*%S3n)
  q14 = -(1/tau)*t(X)%*%S3n
  D = list(NULL)
  for(i in 1:length(K)){
    D[[i]] = alpha[i]*K[[i]]
  }
  q22 = -(1/tau)*t(Reduce("cbind",N))%*%diag(S1n)%*%Reduce("cbind",N)-bdiag(D)
  q23 = list(NULL)
  q24 = list(NULL)
  q32 = list(NULL)
  q42 = list(NULL)
  for(i in 1:length(N)){
    q23[[i]] = -(1/(tau^2))*(t(N[[i]])%*%diag(S1n)%*%(y-mu)-Delta*t(N[[i]])%*%S3n) - alpha[i]*K[[i]]%*%theta[[i]]
    q24[[i]] = -(1/tau)*t(N[[i]])%*%S3n
    q32[[i]] = t(q23[[i]])
    q42[[i]] = t(q24[[i]])
  }
  q33 = (1/(tau^3))*(tau/2-t((y-mu)^2)%*%S1n+2*Delta*t((y-mu))%*%S3n-(Delta^2)*sum(S2n))
  q34 = (1/(tau^2))*(Delta*sum(S2n)-t((y-mu))%*%S3n)
  q44 = -(1/tau)*sum(S2n)
  Q2 = rbind(cbind(q11,Reduce("cbind",q12),q13,q14),
             cbind(Reduce("rbind",q21),q22,Reduce("rbind",q23),
                   Reduce("rbind",q24)),
             cbind(t(q13),Reduce("cbind",q32),q33,q34),
             cbind(t(q14),Reduce("cbind",q42),t(q34),q44))
  if(type == "Ghat"){
    DX = cbind(X,Reduce("cbind",N),rep(0,n),rep(0,n))
    auxq2 = list(NULL)
    for(i in 1:length(N)){
      auxq2[[i]] = (1/tau)*t(N[[i]])%*%diag(S1n)
    }
    Q2y = rbind((1/tau)*t(X)%*%diag(S1n), Reduce("rbind",auxq2),
                (1/tau)*S3n,(1/tau)*(t(y-mu)%*%diag(S1n)-Delta*S3n))
    GH = DX%*%ginv(-as.matrix(Q2))%*%Q2y
    dim(GH)
    hij = diag(GH)
    return(hij)
  }
  if(type == "Gcooks"){
    p = ncol(X)
    Q1 = matrix(0,n,p+ncol(N[[1]])*length(N)+2)
    for(i in 1:n){
      somabeta = 0
      somatheta = list(NULL)
      for(l in 1:length(N)){
        somatheta[[l]] = rep(0,ncol(N[[l]]))
      }
      somatau = 0
      somaDelta = 0
      for(j in 1:n){
        if(j!=i){
          somabeta = somabeta -(1/tau)*(t(X[j,])*Delta*S3n[j]-t(X[j,])*(y[j]-mu[j])*S1n[j])
          for(k in 1:length(N)){
            somatheta[[k]] = somatheta[[k]] -(1/tau)*((N[[k]][j,])*Delta*S3n[j]-(N[[k]][j,])*(y[j]-mu[j])*S1n[j]) - .5*(alpha[[k]]*K[[k]]%*%(theta[[k]]))
          }
          somatau = somatau - 1/(2*tau) + (1/(2*tau^2))*((y[j]-mu[j])^2*S1n[j]-2*(y[j]-mu[j])*Delta*S3n[j]+Delta^2*S2n[j])
          somaDelta = somaDelta - (1/tau)*(Delta*S2n[j]-(y[j]-mu[j])*S3n[j])
        }
      }
      Q1[i,] = c(somabeta, Reduce("c",somatheta), somatau, somaDelta)
    }
    GCD = 0
    for(i in 1:n){
      GCD[i] = t(Q1[i,])%*%ginv(-as.matrix(Q2))%*%Q1[i,]
    }
    return(GCD)
  }
}


LI.ASMCSN = function(obj, type = c("CW","SP,RP","SKP","CP"), xc = NULL){
  X = obj$comp$X
  N = obj$comp$N
  y = obj$y
  K = obj$comp$K
  n = nrow(X)

  alpha = obj$coefficients$alpha
  mu = obj$fitted.values
  beta = obj$coefficients$beta
  theta = obj$coefficients$theta
  sigma2 = obj$coefficients$sigma2
  gama = obj$coefficients$gamma
  nu = obj$coefficients$nu
  g13 = nthroot(gama,3)
  g23 = g13^2
  sigma = sqrt(sigma2)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  del = lambda/sqrt(1+lambda^2)
  Delta = sigma*del/sqrt(1-b^2*del^2)
  tau = sigma2*(1-del^2)/(1-b^2*del^2)
  S1n = obj$expectations$S1n
  S2n = obj$expectations$S2n
  S3n = obj$expectations$S3n
  q11 = -(1/tau)*t(X)%*%diag(S1n)%*%X
  q12 = list(NULL)
  q21 = list(NULL)
  for(i in 1:length(N)){
    q12[[i]] = -(1/tau)*t(X)%*%diag(S1n)%*%N[[i]]
    q21[[i]] = t(q12[[i]])
  }
  q13 = -(1/(tau^2))*(t(X)%*%diag(S1n)%*%(y-mu)-Delta*t(X)%*%S3n)
  q14 = -(1/tau)*t(X)%*%S3n
  D = list(NULL)
  for(i in 1:length(K)){
    D[[i]] = alpha[i]*K[[i]]
  }
  q22 = -(1/tau)*t(Reduce("cbind",N))%*%diag(S1n)%*%Reduce("cbind",N)-bdiag(D)
  q23 = list(NULL)
  q24 = list(NULL)
  q32 = list(NULL)
  q42 = list(NULL)
  for(i in 1:length(N)){
    q23[[i]] = -(1/(tau^2))*(t(N[[i]])%*%diag(S1n)%*%(y-mu)-Delta*t(N[[i]])%*%S3n) - alpha[i]*K[[i]]%*%theta[[i]]
    q24[[i]] = -(1/tau)*t(N[[i]])%*%S3n
    q32[[i]] = t(q23[[i]])
    q42[[i]] = t(q24[[i]])
  }
  q33 = (1/(tau^3))*(tau/2-t((y-mu)^2)%*%S1n+2*Delta*t((y-mu))%*%S3n-(Delta^2)*sum(S2n))
  q34 = (1/(tau^2))*(Delta*sum(S2n)-t((y-mu))%*%S3n)
  q44 = -(1/tau)*sum(S2n)
  Q2 = rbind(cbind(q11,Reduce("cbind",q12),q13,q14),
             cbind(Reduce("rbind",q21),q22,Reduce("rbind",q23),
                   Reduce("rbind",q24)),
             cbind(t(q13),Reduce("cbind",q32),q33,q34),
             cbind(t(q14),Reduce("cbind",q42),t(q34),q44))
  A1 = tau/sigma2
  A2 = -del*sqrt(1-b^2*del^2)/(2*sigma2^(3/2)*(1-del^2))
  A3 = (sigma/sqrt(1-b^2*del^2))*(1+del^2*b^2/(1-b^2*del^2))
  A4 = (2*del*sigma2/(1-b^2*del^2))*((1-del^2)*b^2/(1-b^2*del^2)-1)
  A5 = (tau*A3-Delta*A4)/(tau^2)
  A6 = (2*tau*Delta*A3-Delta^2*A4)/(tau^2)
  A7 = Delta^2*tau/(sigma2^3) - A1*Delta^2/(sigma2^2)
  if(type == "CW"){
    Db = (1/tau)*t(X)%*%(diag((y-mu)*S1n)-diag(S3n)*Delta)
    Dt = list(NULL)
    for(i in 1:length(N)){
      Dt[[i]] = (1/tau)*t(N[[i]])%*%(diag((y-mu)%*%S1n)-diag(S3n)*Delta)
    }
    Ds2 = -1/(2*sigma2)+(1/(2*tau*sigma2))*(y-mu)^2*S1n+A2*(y-mu)*S3n-.5*A7*S2n
    Dg = -(A4/(2*tau))+(A4/(2*tau^2))*(y-mu)^2*S1n+A5*(y-mu)*S3n-.5*A6*S2n
    DDelta = rbind(Db,Reduce("rbind",Dt),t(Ds2),t(Dg))
    qpp = t(DDelta)%*%ginv(-as.matrix(Q2))%*%DDelta
    auxe = eigen(qpp)
    LI = auxe$vectors[,1]
    return(LI)
  }
  if(type == "SP"){
    Db = (1/tau)*t(X)%*%(diag((y-mu)%*%S1n)-.5*diag(S3n)*Delta)
    Dt = list(NULL)
    for(i in 1:length(N)){
      Dt[[i]] = (1/tau)*t(N[[i]])%*%(diag((y-mu)%*%S1n)-.5*diag(S3n)*Delta)
    }
    Ds2 = (A1/(2*tau^2))*(y-mu)^2*S1n + .5*A2*(y-mu)*S3n
    Dg = (A4/(2*tau^2))*(y-mu)^2*S1n + .5*A5*(y-mu)*S3n
    DDelta = rbind(Db,Reduce("rbind",Dt),t(Ds2),t(Dg))
    qpp = t(DDelta)%*%ginv(-as.matrix(Q2))%*%DDelta
    auxe = eigen(qpp)
    LI = auxe$vectors[,1]
    return(LI)
  }
  if(type == "RP"){
    sy = sd(y)
    Db = (sy/tau)*t(X)*S1n
    Dt = list(NULL)
    for(i in 1:length(N)){
      Dt[[i]] = (sy/tau)*t(N[[i]])*S1n
    }
    Ds2 = -(sy/(tau*sigma2))*(y-mu)*S1n+A2*sy*S3n
    Dg = -(A4/(tau^2))*sy*(y-mu)*S1n+A5*sy*S3n
    DDelta = rbind(Db,Reduce("rbind",Dt),t(Ds2),t(Dg))
    qpp = t(DDelta)%*%ginv(-as.matrix(Q2))%*%DDelta
    auxe = eigen(qpp)
    LI = auxe$vectors[,1]
    return(LI)
  }
  if(type == "SKP"){
    Db = del*t(X)*diag(as.matrix(A5*S3n-A4*(y-mu)*S1n/(tau^2)))
    Dt = list(NULL)
    for(i in 1:length(N)){
      Dt[[i]] = del*t(N[[i]])*diag(as.matrix(A5*S3n-A4*(y-mu)*S1n/(tau^2)))
    }
    A2d = -(1/(2*sigma2^(3/2)*(1-del^2)^2))*(sqrt(1-b^2*del^2 + del^2*b^2/sqrt(1-b^2*del^2)))
    A4d = (2*sigma2/((1-b^2*del^2)^2))*((b^2/(1-b^2*del^2))*((1-del^2)*(1-b^2*del^2)+4*del^2*b^2*(1-del^2))-(1+b^2*del^2))
    A3d = (sqrt(sigma2)*del*b^2/((1-b^2*del^2)^(3/2)))*(1+(1/(1-b^2*del^2))*(2*(1-b^2*(del^2))+3*del^2*b^2))
    A5d = (1/(tau^3))*(tau*(tau*A3d-Delta*A4d) - 2*A4d*(tau*A3-Delta*A4))
    A6d = (1/(tau^3))*(tau*(2*A3^2*tau+2*Delta*tau*A3d-Delta^2*A4d)-2*A4*(2*Delta*tau*A3-Delta^2*A4))
    Ds2 = A4*del*(y-mu)^2*S1n/(2*sigma2*tau^2)-del*A2d*(y-mu)*S3n
    Dg = (del/(2*tau^2))*(A4d*tau-A4^2) + (del/(2*tau^4))*(A4d*tau^2-2*A4^2*tau)*(y-mu)^2*S1n - del*A5d*(y-mu)*S3n + del*A6d*S2n/2
    DDelta = rbind(Db,Reduce("rbind",Dt),t(Ds2),t(Dg))
    qpp = t(DDelta)%*%ginv(-as.matrix(Q2))%*%DDelta
    auxe = eigen(qpp)
    LI = auxe$vectors[,1]
    return(LI)
  }
  if(type == "CP"){
    if(is.null(xc)){
      stop("Empty xc")
    }
    er = rep(0,p)
    pe = 1
    er[xc] = 1
    sr = sd(X[,xc])
    Db = (1/tau)*((sr*er%*%t(y-mu)-sr*beta%*%er%*%t(X))%*%diag(S1n)-sr*Delta*er%*%t(S3n))
    Dt = list(NULL)
    for(i in 1:length(N)){
      Dt[[i]] = -(sr/tau)*as.numeric(er%*%beta)*t(N[[i]])%*%diag(S1n)
    }
    Ds2 = (1/(tau*sigma2))*sr*as.numeric(er%*%beta)*(y-mu)*S1n- sr*as.numeric(er%*%beta)*A2*S3n
    Dg = (1/(tau^2))*A4*sr*as.numeric(er%*%beta)*(y-mu)*S1n-sr*as.numeric(er%*%beta)*A5*S3n
    DDelta = rbind(Db,Reduce("rbind",Dt),t(Ds2),t(Dg))
    qpp = t(DDelta)%*%ginv(-as.matrix(Q2))%*%DDelta
    auxe = eigen(qpp)
    LI = auxe$vectors[,1]
    return(LI)
  }
}


