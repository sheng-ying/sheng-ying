#========================================================================================#

# Subject-level data: 
## X1: binary predictor with Pr(X=0)=Pr(X=1)=0.5
## X_2,...X_p: independently follows the standard normal 
## logit Pr(Y=1|X) = -1+0.3X_1+0.2X_2+0.2X_3+0.4X_4+0.4X_5

# The form of the external aggregated information:
## the probabilities of Y=1 given in the subgroups {X_1=1,X_2<=0} and {X_1=1,X_2>0}
## the means of X_1 and X_2

# The density ratio model: 
## Z=(1,X_1,X_2,X_p)
## alpha = (0,0,0,0)

### info is the aggregated information without uncertainty
### info.star is the aggregated information with uncertainty 


#========================================================================================#


source("PELfun 1.r")
source("PELfun 2.r")

set.seed(10)
n <- 200
p <- 20
b0 <- c(-1,0.3,0.2,0.2,0.4,0.4,rep(0,p-5))
kappa <- n/400
setA <- c(2:6) 
setB <- c(1,2,3,p+1)
index <- c(2,3)  
info <- c(0.31,0.38,0.5,0)  
info.star <- c(0.31,0.39,0.48,0.11)
x <- cbind(rep(1,n), rbinom(n,1,0.5), rmvnorm(n,mean = rep(0,p-1),sigma = diag(p-1)))
prob <- apply(x,1,function(u){1-1/(1+exp(b0%*%u))})
y <- sapply(prob,function(u){rbinom(1,1,u)})
g1 <- which(x[,3] <= 0 & x[,2] == 1) 
g2 <- which(x[,3] > 0 & x[,2] == 1)
v <- 1.25
lambda <- lam.star <- tau <- 0.1


# penalized MLE with the elastic net penalty
b.ini <- enlr(x,y)


# the proposed PEL estimator and the corresponding BIC
obj.pel <- pel(x,y,info,lambda,v,b.ini,index)
b.pel <- obj.pel$b.new
xi.pel <- obj.pel$xi.new
bic.pel <- bic(x,y,b.pel,xi.pel,info,index)


# the proposed PEL* estimator and the corresponding BIC
obj.star <- pel(x,y,info.star,lam.star,v,b.ini,index)
b.star <- obj.star$b.new
xi.star <- obj.star$xi.new
bic.star <- bic(x,y,b.star,xi.star,info.star,index)


# the proposed PEL_alpha estimator and the corresponding BIC
obj.ex <- pel.ex(x,y,info,tau,v,setB,b.ini,index)
b.ex <- obj.ex$b.new
a.ex <- obj.ex$a.new
eta <- obj.ex$eta.new
bic.ex <- bicex(x,y,a.ex,b.ex,eta,info,setB,index)


# standard error estimates for nonzero parameters 
cov.obj <- get.cov (kappa, setA, x, b.pel, b.star, info, info.star, index, b.ini, lambda1=lambda, lambda2=lam.star, v1=v, v2=v)
cov.ex <- get.cov.ex (setA, setB, x, b.ini, b.ex, a.ex, info, index, tau, v)
see.pel <- cov.obj$se.pel
see.star <- cov.obj$se.star
see.ex <- cov.ex
  
  