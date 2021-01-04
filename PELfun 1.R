library(Matrix)
library(foreach)
library(mvtnorm)
library(glmnet)


#==============================================================#
# the proposed PEL estimator based on homogeneity assumption
# b: paramter of interest
# xi: Lagrange multipliers 
#==============================================================#


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



# given b, llxi/sxi: negative constraint log full likelihood/score w.r.t. xi

llxi <- function(xi,x,b,info,index){
  prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
  psi <- rbind( prob-info[1], prob-info[2], x[,index[1]]-info[3], x[,index[2]]-info[4])
  psi[1,-g1]<-0; psi[2,-g2]<-0
  xi.psi <- xi%*%psi
  xi.psi[which(xi.psi<(1/n-1))] <- (1/n-1)
  sum(sapply(1+xi.psi,function(u){log(u)}))
}


sxi <- function(xi,x,b,info,index){
  prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
  psi <- rbind( prob-info[1], prob-info[2], x[,index[1]]-info[3], x[,index[2]]-info[4])
  psi[1,-g1]<-0; psi[2,-g2]<-0
  (1/(1+xi%*%psi))%*%t(psi)
}



# given xi, llb/sb: negative constraint log full likelihood/score w.r.t. b 

llb <- function(b,x,y,xi,info,index){
  prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
  prob[which(prob==0)] <- 1/n
  lc <- sum(y*sapply(prob,function(u){log(u)})+(1-y)*sapply(1-prob,function(u){log(u)}))
  lm <- llxi(xi,x,b,info,index)
  -lc + lm
}


sb <- function(b,x,y,xi,info,index){
  prob <- apply(x,1,function(u){1-1/(1+exp(b%*%u))})
  psi <- rbind( prob-info[1], prob-info[2], x[,index[1]]-info[3], x[,index[2]]-info[4])
  psi[1,-g1]<-0; psi[2,-g2]<-0
  d1 <- x*(prob*(1-prob)); d1[-g1,]<-0
  d2 <- x*(prob*(1-prob)); d2[-g2,]<-0
  d <- xi%*%rbind( (1/(1+xi%*%psi))%*%d1, (1/(1+xi%*%psi))%*%d2, matrix(0, nrow=length(index), ncol=ncol(x)) )
  - t(x)%*%(y-prob) + t(d) 
}



# local quadratic approximation (LQA) for alasso

pen <- function(b,b.ini,b.old,lambda,v){
  b.ini[which(b.ini==0)] <- 0.001
  b.old[which(b.old==0)] <- 0.001
  lambda*b[-1]%*%diag(1/abs(2*b.old[-1]*abs(b.ini[-1])^v))%*%b[-1] 
}


dpen <- function(b,b.ini,b.old,lambda,v){
  b.ini[which(b.ini==0)] <- 0.001
  b.old[which(b.old==0)] <- 0.001
  lambda*diag(1/abs(b.old*abs(b.ini)^v))%*%c(0,b[-1]) 
}



# BIC with constant = max(log(log(p)),1)

bic <- function(x,y,b,xi,info,index){
  2*llb(b,x,y,xi,info,index) + sum(b[-1]!=0)*log(n)*max(log(log(p)),1)
}


# iterative algorithm for PEL under homogeneity assumption

pel <- function (x,y,info,lambda,v,b.ini,index){
  K <- length(info)
  xi.old <- xi.ini <- rep(0,K)
  b.old <-  b.ini
  k <- 1
  while (k <= 5){
    
    # given xi, update b
    fb <- function(b){ llb(b,x,y,xi.old,info,index) + pen(b,b.ini,b.old,lambda,v) }
    gb <- function(b){ sb(b,x,y,xi.old,info,index) + dpen(b,b.ini,b.old,lambda,v) }
    b.new <- optim(b.old,fb,gb,method = "L-BFGS-B",control = list(maxit=200))$par
    b.new[abs(b.new) <= 0.001] <- 0
    
    # given b, update xi
    fxi <- function(xi){llxi(xi,x,b.new,info,index)}   
    gxi <- function(xi){sxi(xi,x,b.new,info,index)}
    xi.new <- optim(xi.old,fxi,method = "L-BFGS-B",lower=rep(-1/2,K),upper=rep(1/2,K),control = list(maxit=200,fnscale=-1))$par 
    
    if (max(abs(b.new - b.old)) <= 0.00001){break}
    k <- k+1
    b.old <- b.new
    xi.old <- xi.new
  }
  list(b.new = b.new, xi.new = xi.new)
}




