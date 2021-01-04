library(Matrix)
library(foreach)
library(mvtnorm)
library(glmnet)
library(BB)

#=============================================================================#
# the extended PEL_alpha estimator that accounts for population heterogeneity 
# b: paramter of interest
# a: parameter in the density ratio model 
# eta: Lagrange multipliers 
#=============================================================================#


# elastic net MLE 

enlr <- function (x,y){
  alpha <- c(0:10)/10
  cv.value <- rep(0,length(alpha))
  for (i in 1:length(alpha)){
    cv.obj <- cv.glmnet(x[,-1], y, nfolds=10, family="binomial", alpha=alpha[i])
    cv.value[i] <- min(cv.obj$cvm)
  }
  opt <- which.min(cv.value)
  cv.opt <- cv.glmnet(x[,-1], y, nfolds=10, family="binomial", alpha=alpha[opt])
  obj <- glmnet(x[,-1], y, family="binomial", lambda = cv.opt$lambda.min, alpha=alpha[opt])
  c(obj$a0, as.vector(obj$beta))
}



# given b, obtain initial value for a based on estimation equations

get.a <- function (x,b,info,setB,index){
  fa <- function(a){ 
    xa <- x[,setB]
    psi0 <- t(exp(xa%*%a)-1)
    prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
    psi1 <- exp(xa%*%a)*(prob-info[1]); psi1[-g1]<-0
    psi2 <- exp(xa%*%a)*(prob-info[2]); psi2[-g2]<-0
    psi3 <- exp(xa%*%a)*(x[,index[1]]-info[3])
    psi4 <- exp(xa%*%a)*(x[,index[2]]-info[4])
    mean <- c( mean(psi0), mean(psi1), mean(psi2), mean(psi3), mean(psi4))
    sum (mean^2)
  }
  optim(rep(0,length(setB)),fa,control = list(maxit=100))$par
}


# given a and b, lleta/seta: constraint log full likelihood/score w.r.t. eta (negative)

lleta <- function(eta,x,a,b,info,setB,index){
  xa <- x[,setB]
  psi0 <- t(exp(xa%*%a)-1)
  prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
  psi.ex <- rbind(psi0, as.vector(exp(xa%*%a)*(prob-info[1])), as.vector(exp(xa%*%a)*(prob-info[2])),
                  as.vector(exp(xa%*%a)*(x[,index[1]]-info[3])), as.vector(exp(xa%*%a)*(x[,index[2]]-info[4])) )
  psi.ex[2,-g1]<-0; psi.ex[3,-g2]<-0
  eta.psi <- eta%*%psi.ex
  eta.psi[which(eta.psi<=(1/n-1))] <- (1/n-1)
  eta.psi[which(eta.psi>n)] <- n
  sum(sapply(1+eta.psi,function(u){log(u)}))
}


seta <- function(eta,x,a,b,info,setB,index){
  xa <- x[,setB]
  psi0 <- t(exp(xa%*%a)-1)
  prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
  psi.ex <- rbind(psi0, as.vector(exp(xa%*%a)*(prob-info[1])), as.vector(exp(xa%*%a)*(prob-info[2])),
                  as.vector(exp(xa%*%a)*(x[,index[1]]-info[3])), as.vector(exp(xa%*%a)*(x[,index[2]]-info[4])) )
  psi.ex[2,-g1]<-0; psi.ex[3,-g2]<-0
  (1/(1+eta%*%psi.ex))%*%t(psi.ex)
}


# given eta and a, llb.ex/sb.ex: constraint log full likelihood/score w.r.t. b (negative)

llb.ex <- function(a,b,x,y,eta,info,setB,index){
  prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
  prob[which(prob==0)]<-1/n
  prob[which(prob==1)]<-1-1/n
  lc <- sum(y*sapply(prob,function(u){log(u)})+(1-y)*sapply(1-prob,function(u){log(u)}))
  lm <- lleta(eta,x,a,b,info,setB,index)
  -lc + lm
}


sb.ex <- function(a,b,x,y,eta,info,setB,index){
  xa <- x[,setB]
  psi0 <- t(exp(xa%*%a)-1)
  prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
  psi.ex <- rbind(psi0, as.vector(exp(xa%*%a)*(prob-info[1])), as.vector(exp(xa%*%a)*(prob-info[2])),
                  as.vector(exp(xa%*%a)*(x[,index[1]]-info[3])), as.vector(exp(xa%*%a)*(x[,index[2]]-info[4])) )
  psi.ex[2,-g1]<-0;psi.ex[3,-g2]<-0
  d1.ex <- d4.ex <- d5.ex <- matrix(0,ncol = (p+1),nrow = n)
  d2.ex <- as.vector(exp(xa%*%a))*x*(prob*(1-prob))
  d3.ex <- as.vector(exp(xa%*%a))*x*(prob*(1-prob))
  d2.ex[-g1,]<-0; d3.ex[-g2,]<-0
  d.ex <- eta%*%rbind( (1/(1+eta%*%psi.ex))%*%d1.ex, (1/(1+eta%*%psi.ex))%*%d2.ex, (1/(1+eta%*%psi.ex))%*%d3.ex,
                       (1/(1+eta%*%psi.ex))%*%d4.ex, (1/(1+eta%*%psi.ex))%*%d5.ex)
  - t(x)%*%(y-prob) + t(d.ex) 
}


# local quadratic approximation (LQA) for alasso

pen.ex <- function(b,b.ini,b.old,tau,v){
  b.ini[which(b.ini==0)] <- 0.001
  b.old[which(b.old==0)] <- 0.001
  tau*b[-1]%*%diag(1/abs(2*b.old[-1]*abs(b.ini[-1])^v))%*%b[-1]
}


dpen.ex <- function(b,b.ini,b.old,tau,v){
  b.ini[which(b.ini==0)] <- 0.001
  b.old[which(b.old==0)] <- 0.001
  tau*diag(1/abs(b.old*abs(b.ini)^v))%*%c(0,b[-1]) 
}


# BIC with constant = max(log(log(p)),1)

bicex <- function(x,y,a,b,eta,info,setB,index){
  2*llb.ex(a,b,x,y,eta,info,setB,index) + sum(b[-1]!=0)*log(n)*max(log(log(p)),1)
}


# standard error estimation 

get.cov.ex <- function(setA, setB, x, b.ini, b.ex, a.ex, info, index, tau, v){
  
  x.true <- x[,setA]
  xa <- x[,setB]
  b.ex <- b.ex[setA]
  b1.ex <- sapply(b.ex,function(u){max(u,0.1)})
  b.ini <- sapply(b.ini,function(u){max(u,0.1)})
  weight <- abs(b.ini[setA])^(-v)
  dpen.ex <- tau*diag(weight/b1.ex)/n
  
  prob3 <- apply(x.true,1,function(u){1-1/(1+exp(b.ex%*%u))})
  d5 <- t(x.true)%*%diag(prob3*(1-prob3))%*%x.true
  psi0 <- t(exp(xa%*%a.ex)-1)
  psi.ex <- rbind(psi0, as.vector(exp(xa%*%a.ex)*(prob3-info[1])), as.vector(exp(xa%*%a.ex)*(prob3-info[2])),
                  as.vector(exp(xa%*%a.ex)*(x[,index[1]]-info[3])), as.vector(exp(xa%*%a.ex)*(x[,index[2]]-info[4])) ) 
  psi.ex[2,-g1]<-0; psi.ex[3,-g2]<-0
  cov3.psi <- psi.ex%*%t(psi.ex) 
  d1.ex <- cbind( matrix(0,ncol = (length(setA)),nrow = n),  xa*as.vector(exp(xa%*%a.ex)))
  d2.ex <- cbind( as.vector(exp(xa%*%a.ex))*x.true*(prob3*(1-prob3)), xa*as.vector(exp(xa%*%a.ex)*(prob3-info[1])) ); 
  d3.ex <- cbind( as.vector(exp(xa%*%a.ex))*x.true*(prob3*(1-prob3)), xa*as.vector(exp(xa%*%a.ex)*(prob3-info[2])) ); 
  d2.ex[-g1,]<-0; d3.ex[-g2,]<-0
  d4.ex <- cbind( matrix(0,ncol = length(setA), nrow = n), xa*as.vector(exp(xa%*%a.ex))*x[,index[1]] )
  d5.ex <- cbind( matrix(0,ncol = length(setA), nrow = n), xa*as.vector(exp(xa%*%a.ex))*x[,index[2]] ) 
  d.ex <- rbind(colSums(d1.ex),colSums(d2.ex),colSums(d3.ex),colSums(d4.ex),colSums(d5.ex))
  cov1.ex <- as.matrix(bdiag(d5,matrix(0,ncol=length(setB),nrow=length(setB))))/n + t(d.ex)%*%solve(cov3.psi)%*%d.ex/n
  cov1.ex <- cov1.ex[1:length(setA),1:length(setA)]
  cov.ex <- solve(cov1.ex+dpen.ex)%*%cov1.ex%*%solve(cov1.ex+dpen.ex)
  se.ex <- sqrt(diag(solve(cov.ex))/n)[1:length(setA)]
  return(se.ex)
}



# iterative algorithm for exteneded PEL_alpha 

pel.ex <- function (x,y,info,tau,v,setB,b.ini,index){
  K <- length(info)
  eta.old <- eta.ini <- rep(0,K+1)
  b.old <- b.ini
  a.old <- a.ini <- get.a(x,b.ini,info,setB,index)
  k <- 1
  while (k <= 5){
    
    # given eta and a, update b 
    fb <- function(b){ llb.ex(a.old,b,x,y,eta.old,info,setB,index) + pen.ex(b,b.ini,b.old,tau,v) }
    gb <- function(b){ sb.ex(a.old,b,x,y,eta.old,info,setB,index) + dpen.ex(b,b.ini,b.old,tau,v) }
    b.new <- optim(b.old,fb,gb,method = "L-BFGS-B",control = list(maxit=100))$par
    b.new[which(abs(b.new) < 0.001)] <- 0
    
    # given b and eta, update a
    fa <- function(a){ lleta(eta.old,x,a,b.new,info,setB,index) }   
    a.new <- optim(a.old,fa,method = "L-BFGS-B",control = list(maxit=100))$par
    
    # given a and b, update eta
    feta <- function(eta){ lleta(eta,x,a.new,b.new,info,setB,index) }   
    geta <- function(eta){ seta(eta,x,a.new,b.new,info,setB,index) }
    eta.new <- optim(eta.old,feta,geta,method = "L-BFGS-B",lower=rep(-1/2,(K+1)),upper=rep(1/2,(K+1)),control = list(maxit=100,fnscale=-1))$par
    
    if (max(abs(b.new- b.old)) <= 0.00001){break}
    k <- k+1
    if (llb.ex(a.new,b.new,x,y,eta.new,info,setB,index)< llb.ex(a.old,b.old,x,y,eta.old,info,setB,index)){
      b.old <- b.new
      a.old <- a.new
      eta.old <- eta.new
    }
  }
  list(b.new=b.new, a.new=a.new, eta.new=eta.new)
}


