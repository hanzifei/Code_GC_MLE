library(geoR); library(MASS);library(mvtnorm)
library(Rcpp);library(RcppArmadillo); library(boot)
library(MASS); library(snowfall); library(spatstat)
library(VGAM); library(gcKrig); library(BB)
library(numDeriv); library(pscl); library(gcmr)

###################################################################################################
############# LANSING WOODS Example  ##########
###################################################################################################
data(lansing)
boundx1 <- rep(seq(0, 15/16, by = 1/16), 16)
boundx2 <- rep(seq(0, 15/16, by = 1/16), 16)+1/16
boundy1 <- rep(seq(0, 15/16, by = 1/16), each = 16)
boundy2 <- rep(seq(0, 15/16, by = 1/16), each = 16)+1/16
hickory <- c()
for (i in 1:256){
  hickory[i] = length(which(lansing$x > boundx1[i] & lansing$x <= boundx2[i]
                            & lansing$y >= boundy1[i] & lansing$y < boundy2[i]
                            &    lansing$marks=='hickory'))
}

maple <- c()
for (i in 1:256){
  maple[i] = length(which(lansing$x > boundx1[i] & lansing$x <= boundx2[i]
                          & lansing$y >= boundy1[i] & lansing$y < boundy2[i]
                          &    lansing$marks=='maple'))
}

redoak <- c()
for (i in 1:256){
  redoak[i] = length(which(lansing$x > boundx1[i] & lansing$x <= boundx2[i]
                           & lansing$y >= boundy1[i] & lansing$y < boundy2[i]
                           &    lansing$marks=='redoak'))
}

whiteoak <- c()
for (i in 1:256){
  whiteoak[i] = length(which(lansing$x > boundx1[i] & lansing$x <= boundx2[i]
                             & lansing$y >= boundy1[i] & lansing$y < boundy2[i]
                             &    lansing$marks=='whiteoak'))
}
oaks = redoak+whiteoak
rm(i)

###### Make Hickory and Maple as Geostatistical Count Data 
wood.xloc <- rep(seq(1/32, 31/32, by = 2/32), 16)
wood.yloc <- rep(seq(1/32, 31/32, by = 2/32), each = 16)
distwood <- as.matrix(dist(cbind(wood.xloc,wood.yloc),method='euclidean'))
hickorymat <- cbind(wood.xloc, wood.yloc, hickory)
maplemat <- cbind(wood.xloc, wood.yloc, maple)
twotree <- cbind(wood.xloc, wood.yloc, hickory, maple,redoak, whiteoak, oaks)
hickx <- twotree[,3]; mapley <- twotree[,4]

###################################################################################################
################### Analysis of the data using GHK method via package gcmr  ######################
################## To do so, need to modify one function in the package ######################
###################################################################################################
 matern_nugget.cormat <- function (D, alpha = 0.5) 
 {
   ans <- list()
   ans$npar <- 2
   ans$start <- function() {
     tau <- c(0.2, 0.2) 
     names(tau) <- c("range", "nugget")
     attr(tau, "upper") <- c(400,1)
     attr(tau, "lower") <- rep(sqrt(.Machine$double.eps), 2)
     tau
   }
   ans$chol <- function(tau, not.na) {
     S <- (1 - tau[2]) * geoR::matern(D, tau[1], alpha) + tau[2] * diag(NROW(D))
     q <- try(chol(S[not.na, not.na]), silent = TRUE)
     if (inherits(q, "try-error")) 
       NULL
     else q
   }
   class(ans) <- c("matern_nugget.gcmr", "cormat.gcmr")
   ans
 }


### GHK estimates via package gcKrig

estpois00 <- mlegc(y = LansingTrees$maple, x = LansingTrees$hickory, 
                   locs = cbind(LansingTrees$Easting, LansingTrees$Northing), 
                   marginal = poisson.gc(link = 'log'), 
                   corr = matern.gc(nugget = TRUE))

###################################################################################################
####################################### DC Algorithm #################################################
###################################################################################################
############# Prep Functions
#### Functions Required in the process of DC. 
#### To avoid overflow
dlmvnormx <- function(x,mean,sigma){
  llik <- dmvnorm(x, mean = mean, sigma = sigma, log = T)
  llik <-  ifelse(llik==-Inf, -.Machine$double.xmax, llik)
  llik <- ifelse(llik== Inf, .Machine$double.xmax, llik)
  llik[is.nan(llik)] == log(.Machine$double.xmin)
  return(sum(llik))
}

#### Log likelihood of: prod(hi(yi|latent, parameters)) with cloned latents
ppois1 <- function(data, tausq, gammas, mu){
  prob = (pnorm((qnorm(ppois(q = data, lambda = mu))-gammas)/sqrt(tausq))-
            pnorm((qnorm(ppois(q = data-1, lambda = mu))-gammas)/sqrt(tausq)))
  prob2 = ifelse(prob==0, .Machine$double.eps, prob)
  sumloglike2 = sum(log(prob2))
  return(sumloglike2)
}

#### Returns the Negative Loglikelihood of the data and augmented data (latent fields)
llik1 <- function(x, mu, data, theta, tausq){
  prob =  (pnorm((qnorm(ppois(q = data, lambda = mu))-x)/sqrt(tausq))-
             pnorm((qnorm(ppois(q = data-1, lambda = mu))-x)/sqrt(tausq)))
  prob = ifelse(prob==0, .Machine$double.eps, prob)
  sumloglike1 = sum(log(prob))
  prob2 = dmvnorm(x = x, mean = rep(0, length(x)), sigma = (1-tausq)*exp(-distwood/theta), log = T)
  return(-(prob2+sumloglike1))
}

#### Fast Rcpp functions to compute: 
#### Generate multivariate normal r.v. matrix with clones (so the means may different )
#### with same correlation parameters 

cppFunction(depends = 'RcppArmadillo',   code = '
            arma::mat matmvnrmCpp (arma::mat mu, arma::mat Sigma){
            int k = mu.n_cols;
            int n = mu.n_rows;
            arma::mat Y = arma::randn(k, n);
            return (mu.t() + Y*arma::chol(Sigma)).t(); 
            }')

#### Compute log-likelihood of mvn matrix, with different mean vector. 

cppFunction(depends = 'RcppArmadillo',   code = '
            double dlmvnmatCpp(arma::mat x, arma::mat mu, arma::mat Sigma ) { 
            const double log2pi = std::log(2.0 * M_PI);
            int n = x.n_cols; int xdim = x.n_rows;
            arma::vec out(n);
            double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
            arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
            double rootisum = arma::sum(log(rooti.diag()));
            for (int i=0; i < n; i++) {
            arma::vec z = rooti * (x.col(i) - mu.col(i));
            out(i) = constants - 0.5 * arma::sum(z%z) + rootisum;
            if(out(i) == 0 ){
            out(i) =  std::numeric_limits<float>::min() ;
            }
            }
            return arma::sum(out);
            }')


#### Gradient of the latent field...... of the Joint log-likelihood
xi <- function(x, mu) {
  out = qnorm(ppois(x, mu))
  out[is.infinite(out)] = 0
  return(out)
}

grad.joint <- function(x, mu, data, theta, tausq){
  # x = gammanew; mu = munew; data = mapley; theta = theta0; tausq = tausq0
  num = dnorm((qnorm(ppois(data-1,mu))-x)/sqrt(tausq)) - dnorm((qnorm(ppois(data,mu))-x)/sqrt(tausq))
  #256*k
  denum = pnorm((qnorm(ppois(data,mu))-x)/sqrt(tausq)) - pnorm((qnorm(ppois(data-1,mu))-x)/sqrt(tausq))
  denum[denum ==0] = .Machine$double.eps
  #256*k
  gradmvn = solve( (1-tausq)*exp(-distwood/theta))%*%(-x) 
  return(gradmvn + num/(denum*sqrt(tausq)))
}

####  Gradient of the beta1, beta2 in log-posterior respectively

betapost1 <- function(beta, data, tausq, gammas, beta2){
  prob = (pnorm((qnorm(ppois(q = data, lambda = exp(beta+beta2*hickx)))-gammas)/sqrt(tausq))-
            pnorm((qnorm(ppois(q = data-1, lambda = exp(beta+beta2*hickx)))-gammas)/sqrt(tausq)))
  prob2 = ifelse(prob==0, .Machine$double.eps, prob)
  sumloglike2 = sum(log(prob2))
  return(sumloglike2)
}

betapost2 <- function(beta, data, tausq, gammas, beta1){
  prob = (pnorm((qnorm(ppois(q = data, lambda = exp(beta1+beta*hickx)))-gammas)/sqrt(tausq))-
            pnorm((qnorm(ppois(q = data-1, lambda = exp(beta1+beta*hickx)))-gammas)/sqrt(tausq)))
  prob2 = ifelse(prob==0, .Machine$double.eps, prob)
  sumloglike2 = sum(log(prob2))
  return(sumloglike2)
}


post.theta <- function(theta, data, gammas, tausq){
  den = dlnorm(theta, mean = -3, sd= 1, log = T)+
    dlmvnormx(t(gammas),mean=rep(0,256),sigma = (1-tausq)*exp(-distwood/theta))
  return(den)
}


grad.logit.tausq <- function(data, gammas, mu, theta, tausq){
  a = dnorm((qnorm(ppois(data-1, mu)) - gammas)/
              sqrt(tausq))*(qnorm(ppois(data-1, mu))-gammas)
  a[is.nan(a)] = 0
  b = dnorm((qnorm(ppois(data, mu)) - gammas)/
              sqrt(tausq))*(qnorm(ppois(data, mu))-gammas)
  b[is.nan(b)] = 0
  c = ((pnorm((qnorm(ppois(q = data, lambda = mu))-gammas)/sqrt(tausq))-
          pnorm((qnorm(ppois(q = data-1, lambda = mu))-gammas)/sqrt(tausq))))
  c[c==0] = .Machine$double.eps
  d = 0.5*sum((a-b)/c)*tausq^{-1.5}
  e = 256*k/(2*(1-tausq)) - sum(diag(t(gammas)%*%solve(exp(-distwood/theta))%*%gammas)/(2*(1-tausq)^2))
  return((d+e)*(tausq)*(1-tausq))
}


#### Function in calculating the variance of proposal distribution of 
#### Gaussian approximation of proposal using Taylor Expansion up to 2.
#### This is used only once for getting size of the jump approximately (for time reason) 

propvarfun <- function(x, mu, data, theta, tausq){
  atmp = ( ((x-xi(data,mu))/sqrt(tausq))*dnorm((xi(data,mu)-x)/sqrt(tausq)) - 
             ((x-xi(data-1,mu))/sqrt(tausq))*dnorm((xi(data-1,mu)-x)/sqrt(tausq)) )/
    (pnorm((xi(data,mu)-x)/sqrt(tausq)) - pnorm((xi(data-1,mu)-x)/sqrt(tausq)))
  
  btmp = (( dnorm((xi(data-1,mu)-x)/sqrt(tausq))  -  dnorm((xi(data,mu)-x)/sqrt(tausq)) )/
            ( pnorm((xi(data,mu)-x)/sqrt(tausq))  -  pnorm((xi(data-1,mu)-x)/sqrt(tausq)) ))^2
  
  c2 = (atmp-btmp)/tausq
  corr1 = solve((1-tausq)*exp(-distwood/theta))
  cov2 = (solve(corr1-diag(c2)))
  return(cov2)
}


#### Initial Values

beta0 = estpois00$MLE[1:2]
theta0 = estpois00$MLE[3]
tausq0 = estpois00$MLE[4]
# beta0 = c(1.11,-0.15);theta0 = 0.08; tausq0 = 0.51
#### Initial Value of gamma0
z01 = qnorm(ppois(mapley-1, lambda = exp(beta0[1]+beta0[2]*hickx)))
z02 = qnorm(ppois(mapley, lambda = exp(beta0[1]+beta0[2]*hickx)))
z0 = (dnorm(z01)-dnorm(z02))/(pnorm(z02)-pnorm(z01))

#### Posterior mode
tmpPois = BB::spg(par = z0, fn = llik1, mu = exp(beta0[1]+beta0[2]*hickx), data = mapley,
              theta = theta0, tausq = tausq0, control = c(maxit = 10000, ftol = 1e-5, gtol = 1e-5,trace = F))

#### Posterior Var if using normal approximation with Taylor Expansion up to 2
propmu = tmpPois$par
propvar = propvarfun(x = tmpPois$par, mu = exp(beta0[1]+beta0[2]*hickx), data = mapley,
                     theta = theta0, tausq = tausq0)

#####################################################################################################################
############# SET UP THE CHAIN
#####################################################################################################################
nit <- 5000; #nit <- 15000
k = 40
theta <- rep(0, nit); tausq <- rep(0, nit)
beta <- matrix(0, nrow = nit, ncol = 2)
gamma <- array(0, dim = c(256, k, nit))
mu <- matrix(0, nrow = nit, ncol = 256)

gamma0 <- t(mvrnorm(k, mu = propmu, Sigma = propvar)) 

# Initialize of all parameters including the latent fields with k copies
theta[1] <- theta0; tausq[1] <- tausq0
beta[1,] <- beta0;  gamma[,,1] <- gamma0
mu[1,] <- exp(beta0[1] + beta0[2]*hickx)

####
#### Start the Chain
#### Two Blocks Langevin Updates

counter1 <- 0; counter2 <- 0

#########################################################################################################################
#########################################################################################################################

h1 = 0.000085; h2 = 0.000003; h3 = 0.09
# Langevin stepsize of beta1, beta2, gamma
tt = 0.00024
gg = 0.00088 

#trans.nod <- seq(500,nit,by = 500)  
# Where Adaptive may Take Place

for (i in 1:(nit-1)){
  
  # For RWP we can use this 
  #  if (i %in% trans.nod){
  #    step.beta = as.vector(cov(beta[(i-499):i,]))*(2.38/sqrt(2))^2
  #    step.theta = sd(log(theta[(i-499):i]))*2.38/sqrt(2)
  #    step.other.sd[2] = sd(logit(tausq[(i-499):i]))*2.38/sqrt(2)
  #  }
  
  #### Calculate Gradient of beta1 and beta2 and
  #### Set Langevin Updates
  
  grad11 = numDeriv::grad(function(x) {betapost1(x, data = mapley, tausq = tausq[i], 
                                       gammas = gamma[,,i], beta2 = beta[i,2])},beta[i,1]) - (beta[i,1])/30
  
  bmunew1 = beta[i,1] + (h1/2)*grad11
  betanew1 <- rnorm(1, mean = bmunew1, sd = sqrt(h1))
  grad12 = numDeriv::grad(function(x) {betapost1(x, data = mapley, tausq = tausq[i], 
                                       gammas = gamma[,,i], beta2 = beta[i,2])},betanew1) - (betanew1)/30
  bmu1back <- betanew1+(h1/2)*grad12
  grad21 =  numDeriv::grad(function(x) {betapost2(x, data = mapley, tausq = tausq[i], 
                                        gammas = gamma[,,i], beta1 = beta[i,1])},beta[i,2]) - (beta[i,2])/30
  
  bmunew2 = beta[i,2] + (h2/2)*grad21
  betanew2 <- rnorm(1, mean = bmunew2, sd = sqrt(h2))
  munew <- exp(betanew1+betanew2*hickx)

  grad22 =  grad(function(x) {betapost2(x, data = mapley, tausq = tausq[i], 
                                        gammas = gamma[,,i], beta1 = beta[i,1])},betanew2) - (betanew2)/30

  bmu2back <- betanew2 + (h2/2)*grad22
  
  logratio1 <- min(0, (dnorm(betanew1, mean = 0, sd = sqrt(30), log = T)+
                         dnorm(betanew2, mean = 0, sd = sqrt(30), log = T)-
                         dnorm(beta[i,1], mean = 0, sd = sqrt(30), log = T)-
                         dnorm(beta[i,2], mean = 0, sd = sqrt(30), log = T)+
                         ppois1(mapley, tausq[i], gamma[,,i], munew)- 
                         ppois1(mapley, tausq[i], gamma[,,i], mu[i, ])+
                         dnorm(beta[i, 1], mean = bmunew1, sd = sqrt(h1), log = T) -
                         dnorm(betanew1, mean = bmu1back, sd = sqrt(h1), log = T) +
                         dnorm(beta[i, 2], mean = bmunew2, sd = sqrt(h2), log = T) -
                         dnorm(betanew2, mean = bmu2back, sd = sqrt(h2), log = T) 
  ))
  
  if (runif(1) <= exp(logratio1)){
    beta[i+1,1] = betanew1; beta[i+1,2] = betanew2; mu[i+1, ] = munew; counter1 = counter1+1
  } else {
    beta[i+1,1] = beta[i,1]; beta[i+1,2] = beta[i,2]
    mu[i+1,] = mu[i,]; counter1 = counter1 
  }
  
  #### Second Block: Update latent fields and the dependence parameters  
  #### For tausq: using Langevin Update to logit(tausq)
  #### For theta: using rw to log(theta)
  #### Tentatively work fine, may be improved anyway.
  
    gradltheta <- grad(function(x) {post.theta(x, data = mapley, gammas = gamma[,,i], 
                                               tausq = tausq[i])},theta[i])*theta[i]
    ltheta <- log(theta[i])
    multhetanew <- ltheta+0.5*tt*gradltheta
    lthetanew <- rnorm(1, mean = multhetanew, sd = sqrt(tt))
    thetanew <- exp(lthetanew)
    
    gradlthetanew <- grad(function(x) {post.theta(x, data = mapley, gammas = gamma[,,i], 
                                               tausq = tausq[i])},thetanew)*thetanew
    multhetaback <- lthetanew+0.5*tt*gradlthetanew
    
  phi <- logit(tausq[i])
  muphinew <- phi+0.5*gg*grad.logit.tausq(mapley, gamma[,,i], mu[i+1, ], theta[i], tausq[i])
  phinew <- rnorm(1, mean = muphinew, sd = sqrt(gg))
  
  tausqnew <- inv.logit(phinew)
  tausqnew <- ifelse(tausqnew > 0.99, 0.99, tausqnew) 
  tausqnew <- ifelse(tausqnew < 0.0001, 0.0001, tausqnew) 
  
  muphiback <- phinew + 0.5*gg*grad.logit.tausq(mapley, gamma[,,i], mu[i+1, ], theta[i], tausqnew)
  grad.mat1 <- grad.joint(x = gamma[,,i], mu = mu[i+1,], data = mapley, theta = theta[i], tausq = tausq[i])
  mu.gammanew = gamma[,,i] + (h3/2)*propvar%*%grad.mat1
  gammanew = matmvnrmCpp(mu = mu.gammanew, Sigma = h3*propvar)
  grad.mat2 <- grad.joint(x = gammanew, mu = mu[i+1,], data = mapley, theta = theta[i], tausq = tausq[i])
  mu.gammaback = gammanew + (h3/2)*propvar%*%grad.mat2
  
  logratio2 <- min(0, (    ppois1(mapley,tausqnew, gammanew, mu[i+1,])- 
                           ppois1(mapley, tausq[i], gamma[,,i], mu[i+1, ])+
                           dlmvnormx(t(gammanew),mean=rep(0,256),sigma = (1-tausqnew)*exp(-distwood/thetanew)) -
                           dlmvnormx(t(gamma[,,i]),mean=rep(0,256),sigma = (1-tausq[i])*exp(-distwood/theta[i])) +
                           dnorm(phi, mean = muphinew, sd = sqrt(gg), log = T) -
                           dnorm(phinew, mean = muphiback, sd = sqrt(gg), log = T) +
                           dnorm(ltheta, mean = multhetanew, sd = sqrt(tt), log = T) -
                           dnorm(lthetanew, mean = multhetaback, sd = sqrt(tt), log = T)+
                           log(tausqnew)+log(1-tausqnew)-log(tausq[i])-log(1-tausq[i])+log(thetanew)-log(theta[i]) + 
                             #The Jacobian
                           dlnorm(theta[i], mean = -2, sd= 1, log = T)-dlnorm(thetanew, mean = -2, sd= 1, log = T) + 
                             #Arbitrary Prior
                           dlmvnmatCpp(gamma[,,i], mu = mu.gammaback, Sigma = h3*propvar)-
                           dlmvnmatCpp(gammanew, mu = mu.gammanew, Sigma = h3*propvar)
  ))
  if (runif(1) <= exp(logratio2)){
    tausq[i+1] = tausqnew; theta[i+1] = thetanew; gamma[,,i+1] = gammanew; counter2 = counter2+1
  } else {
    theta[i+1] = theta[i]; tausq[i+1] = tausq[i]; gamma[,,i+1] = gamma[,,i]; counter2 = counter2
  }
  
  #### Save too many gamma will slow down your computer and eat lots of memory  !!!
  #### If delete some gamma during the process will faster the compuation
  #### 10000 iters will eat 1.1GB for gamma.
  #### Delete one by one will resolve the issue and increase speed.
  #####
}
counter1/nit;  counter2/nit

#sqrt(k*(var(beta[5000:nit,1]))) #0.123
#sqrt(k*var(beta[5000:nit,2])) #0.025
#sqrt(k*var(theta[5000:nit])) # 0.018
#sqrt(k*var(tausq[5000:nit])) #0.101
#plot(gamma[2,1,],type='l')
#acf(gamma[1,3,5000:10000])
#mean(beta[5000:nit,1]) # 1.123
#mean(beta[5000:nit,2]) #-0.150
#mean(theta[5000:nit]) #0.083
#mean(tausq[5000:nit]) #0.494

  #####################################################################################################################
 ############# Method of the Nikoloulopoulos (2011, 2015) using Genz and Bretz via mvtnorm
 #####################################################################################################################
 
 directlik <- function(v)
 {
   #v = c(1.1,-0.15,0.08,0.5)
   beta0 = v[1]; beta1 = v[2]; theta = v[3]; tausq = v[4]
   mu = exp(beta0 + beta1*hickx)
   
   lower = qnorm(ppois(q = mapley-1, lambda = mu, lower.tail = TRUE, log.p = FALSE))
   upper = qnorm(ppois(q = mapley, lambda = mu, lower.tail = TRUE, log.p = FALSE))
   
   R = (1-tausq)*exp(-distwood/theta)+(distwood==0)*tausq
   
   set.seed(1234)
   lik = -log(pmvnorm(lower = lower, upper = upper, mean = rep(0, length(lower)), sigma = R, 
                      algorithm = GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0))[1])
   
   lik = ifelse(lik == Inf, -log(.Machine$double.xmin), lik)
   lik = ifelse(lik == -Inf, log(.Machine$double.xmin), lik)
   return(lik)
 }
 
 gbresult <- try(optim(par = c(1.2,-0.2,0.07,0.5), fn = directlik, gr = NULL, method = c("L-BFGS-B"),
                     lower = c(-3, -3, 0.01, 0.01), upper = c(5, 5, 4, 0.95), hessian = T), silent = T)
 
 
 gbresult.gb <- optim(par = c(1.1,-0.15,0.08,0.5), fn = directlik, method = "L-BFGS-B", 
                      lower = c(-Inf, -Inf, 0, 0), upper = c(Inf, Inf, Inf, 0.9), hessian = T)
 #[1]  1.10011017 -0.14932924  0.08000855  0.49998692
# hess(gbresult.gb)
# sqrt(diag(solve(gbresult.gb$hessian)))
 
 
 
 
################################################################################################################
################################################################################################################
############################ Analyze the Lansing Woods Data with ZIP Marginal ##################################
################################################################################################################
################################################################################################################

 
 #### The ZIP Marginal
 ## using other packages
 dzip2 <- function(x, size, mu)
 {
   lambda=mu+mu/size
   pi=1/(size+1)
   out = dzipois(x=x, lambda=lambda, pstr0 = pi)
   return(out)
 }
 pzip2 <- function(q, size, mu)
 {
   lambda=mu+mu/size
   pi=1/(size+1)
   out = pzipois(q=q, lambda=lambda, pstr0 = pi)
   return(out)
 }
 qzip2 <- function(p, size, mu)
 {
   lambda=mu+mu/size
   pi=1/(size+1)
   out = qzipois(p=p, lambda=lambda, pstr0 = pi)
   return(out)
 }
 
 #### Vectorized Original VERSION
 qzip <- function(p, od, mu){
   qpois(pmax( 0, (p-od/(1+od))/(1-od/(1+od)) ), (1+od)*mu)
 }
 
 pzip <- function(q, od, mu){
   (q >= 0)*od/(1+od)+ppois(q, (1+od)*mu)/(1+od)
 }
 ####
 
 
 #### Specify a zero-inflated Poisson marginal in package gcmr 
 
 zip2.marg <- function (link = "log") 
 {
   ans <- list()
   ans$start <- function(y, x, z, offset) {
     if (!is.null(z)) 
       offset <- list(as.vector(offset$mean), as.vector(offset$precision))
     eps <- sqrt(.Machine$double.eps)
     m <- zeroinfl(y ~ x[,2]|1, offset = offset$mean)
     mu <- fitted(m)
     kappa <- max(10 * eps, mean(((y - mu)^2 - mu)/mu^2))
     lambda <- c(coef(m)[1:NCOL(x)], rep.int(0, NCOL(z)))
     lambda[NCOL(x) + 1] <- ifelse(is.null(z), kappa, log(kappa))
     if (is.null(z)) {
       names(lambda) <- c(dimnames(as.matrix(x))[[2L]], 
                          "dispersion")
       attr(lambda, "lower") <- c(rep(-Inf, NCOL(x)), eps)
     }
     else names(lambda) <- c(paste("mean", dimnames(as.matrix(x))[[2L]], 
                                   sep = "."), paste("dispersion", dimnames(as.matrix(z))[[2L]], 
                                                     sep = "."))
     lambda
   }
   ans$npar <- function(x,z)
   {ifelse(!is.null(z), NCOL(x) + NCOL(z), NCOL(x)+1)}
   ans$dp <- function(y, x, z, offset, lambda) {
     nb <- length(lambda)
     mu <- exp( x %*% lambda[1:NCOL(x)] + offset$mean)
     if(is.null(z))
       size <- 1/lambda[nb]
     else 
       size <- 1/exp(z%*%lambda[(NCOL(x)+1):nb])
     cbind(dzip2(y, size = size, mu = mu), pzip2(y, size = size, mu = mu))
   }
   ans$q <- function(p, x, z, offset, lambda) {
     nb <- length(lambda)
     mu <- exp( x %*% lambda[1:NCOL(x)]+ offset$meaen)
     if(is.null(z))
       size <- 1/lambda[nb]
     else 
       size <- 1/exp(z%*%lambda[(NCOL(x)+1):nb])
     qzip2(p,size = size, mu = mu)
   }
   ans$fitted.val <- function(x, z, offset, lambda) {
     exp( x %*% lambda[1:NCOL(x)]+offset$mean)
   }
   ans$type <- "integer"
   class(ans) <- c("marginal.gcmr")
   ans
 }
 
 matern_nugget.cormat <- function (D, alpha = 0.5) 
 {
   ans <- list()
   ans$npar <- 2
   ans$start <- function() {
     tau <- c(0.2, 0.2) 
     names(tau) <- c("range", "nugget")
     attr(tau, "upper") <- c(400,1)
     attr(tau, "lower") <- rep(sqrt(.Machine$double.eps), 2)
     tau
   }
   ans$chol <- function(tau, not.na) {
     S <- (1 - tau[2]) * geoR::matern(D, tau[1], alpha) + tau[2] * diag(NROW(D))
     q <- try(chol(S[not.na, not.na]), silent = TRUE)
     if (inherits(q, "try-error")) 
       NULL
     else q
   }
   class(ans) <- c("matern_nugget.gcmr", "cormat.gcmr")
   ans
 }
 
 ###### GHK Method
 
 estzip00 <- gcmr(mapley~hickx, marginal = zip2.marg(link='log'), cormat = matern_nugget.cormat(distwood),
                  options = gcmr.options(seed = 1234, nrep = c(200,5000)))
 
 summary(estzip00)
 
 ###################################################################################################
 ####################################### DC Algorithm #################################################
 ###################################################################################################
 ############# Prep Functions
 #### Functions Required in the process of DC. 
 #### To avoid overflow
 dlmvnormx <- function(x,mean,sigma){
   llik <- dmvnorm(x, mean = mean, sigma = sigma, log = T)
   llik <-  ifelse(llik==-Inf, -.Machine$double.xmax, llik)
   llik <- ifelse(llik== Inf, .Machine$double.xmax, llik)
   llik[is.nan(llik)] == log(.Machine$double.xmin)
   return(sum(llik))
 }
 
 pzip1 <- function(data, tausq, gammas, mu, od){
   prob = (pnorm((qnorm(pzip(q = data, od = od, mu = mu))-gammas)/sqrt(tausq))-
             pnorm((qnorm(pzip(q = data-1, od = od, mu = mu))-gammas)/sqrt(tausq)))
   prob2 = ifelse(prob==0, .Machine$double.eps, prob)
   sumloglike2 = sum(log(prob2))
   return(sumloglike2)
 }
 
 llik11 <- function(x, mu, od, data, theta, tausq){
   prob =  (pnorm((qnorm(pzip(q = data, od = od, mu = mu))-x)/sqrt(tausq))-
              pnorm((qnorm(pzip(q = data-1, od = od, mu = mu))-x)/sqrt(tausq)))
   prob = ifelse(prob==0, .Machine$double.eps, prob)
   sumloglike1 = sum(log(prob))
   prob2 = dmvnorm(x = x, mean = rep(0, length(x)), sigma = (1-tausq)*exp(-distwood/theta), log = T)
   return(-(prob2+sumloglike1))
 }
 
 #### Fast Rcpp functions to compute: 
 #### Generate multivariate normal r.v. matrix with clones (so the means may different )
 #### with same correlation parameters 
 cppFunction(depends = 'RcppArmadillo',   code = '
             arma::mat matmvnrmCpp (arma::mat mu, arma::mat Sigma){
             int k = mu.n_cols;
             int n = mu.n_rows;
             arma::mat Y = arma::randn(k, n);
             return (mu.t() + Y*arma::chol(Sigma)).t(); 
             }')

 #### Compute log-likelihood of mvn matrix, with different mean vector. 
 cppFunction(depends = 'RcppArmadillo',   code = '
             double dlmvnmatCpp(arma::mat x, arma::mat mu, arma::mat Sigma ) { 
             const double log2pi = std::log(2.0 * M_PI);
             int n = x.n_cols; int xdim = x.n_rows;
             arma::vec out(n);
             double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
             arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
             double rootisum = arma::sum(log(rooti.diag()));
             for (int i=0; i < n; i++) {
             arma::vec z = rooti * (x.col(i) - mu.col(i));
             out(i) = constants - 0.5 * arma::sum(z%z) + rootisum;
             if(out(i) == 0 ){
             out(i) =  std::numeric_limits<float>::min() ;
             }
             }
             return arma::sum(out);
             }')

 k = 40
 
 xi2 <- function(x, mu, od) {
   out = qnorm(pzip(q = x, od = od, mu = mu))
   out[is.infinite(out)] = 0
   return(out)
 }
 
 grad.joint2 <- function(x, mu, od, data, theta, tausq){
   num = dnorm((qnorm(pzip(q = data-1, od = od, mu = mu))-x)/sqrt(tausq)) - 
     dnorm((qnorm(pzip(q = data, od = od, mu = mu))-x)/sqrt(tausq))
   denum = pnorm((qnorm(pzip(q = data, od = od, mu = mu))-x)/sqrt(tausq)) - 
     pnorm((qnorm(pzip(q = data-1, od = od, mu = mu))-x)/sqrt(tausq))
   denum[denum ==0] = .Machine$double.eps
   gradmvn = solve( (1-tausq)*exp(-distwood/theta))%*%(-x) 
   return(gradmvn + num/(denum*sqrt(tausq)))
 }
 
 ####  Gradient of the beta1, beta2 in log-posterior respectively
 betapost12 <- function(beta, od, data, tausq, gammas, beta2){
   prob = (pnorm((qnorm(pzip(q = data, od = od, mu = exp(beta+beta2*hickx)))-gammas)/sqrt(tausq))-
             pnorm((qnorm(pzip(q = data-1, od = od, mu = exp(beta+beta2*hickx)))-gammas)/sqrt(tausq)))
   prob2 = ifelse(prob==0, .Machine$double.eps, prob)
   sumloglike2 = sum(log(prob2))
   return(sumloglike2)
 }
 
 betapost22 <- function(beta, od, data, tausq, gammas, beta1){
   prob = (pnorm((qnorm(pzip(q = data, od = od, mu = exp(beta1+beta*hickx)))-gammas)/sqrt(tausq))-
             pnorm((qnorm(pzip(q = data-1, od = od, mu = exp(beta1+beta*hickx)))-gammas)/sqrt(tausq)))
   prob2 = ifelse(prob==0, .Machine$double.eps, prob)
   sumloglike2 = sum(log(prob2))
   return(sumloglike2)
 }
 
 odpost <- function(od, beta1, beta2, data, tausq, gammas){
   prob = (pnorm((qnorm(pzip(q = data, od = od, mu = exp(beta1+beta2*hickx)))-gammas)/sqrt(tausq))-
             pnorm((qnorm(pzip(q = data-1, od = od, mu = exp(beta1+beta2*hickx)))-gammas)/sqrt(tausq)))
   prob2 = ifelse(prob==0, .Machine$double.eps, prob)
   sumloglike2 = sum(log(prob2))
   return(sumloglike2)
 }
 
 post.theta <- function(theta, data, gammas, tausq){
   den = dlnorm(theta, meanlog = -3, sdlog = 1, log = T)+
     dlmvnormx(t(gammas),mean=rep(0,256),sigma = (1-tausq)*exp(-distwood/theta))
   return(den)
 }
 
 grad.logit.tausq2 <- function(data, od, gammas, mu, theta, tausq){
   a = dnorm((qnorm(pzip(q = data-1, od = od, mu = mu)) - gammas)/
               sqrt(tausq))*(qnorm(pzip(q = data-1, od = od, mu = mu))-gammas)
   a[is.nan(a)] = 0
   b = dnorm((qnorm(pzip(q = data, od = od, mu = mu)) - gammas)/
               sqrt(tausq))*(qnorm(pzip(q = data, od = od, mu = mu))-gammas)
   b[is.nan(b)] = 0
   c = ((pnorm((qnorm(pzip(q = data, od = od, mu = mu))-gammas)/sqrt(tausq))-
           pnorm((qnorm(pzip(q = data-1, od = od, mu = mu))-gammas)/sqrt(tausq))))
   c[c==0] = .Machine$double.eps
   d = 0.5*sum((a-b)/c)*tausq^{-1.5}
   e = 256*k/(2*(1-tausq)) - sum(diag(t(gammas)%*%solve(exp(-distwood/theta))%*%gammas)/(2*(1-tausq)^2))
   return((d+e)*(tausq)*(1-tausq))
 }
 
 
 #### Function in calculating the variance of proposal distribution of 
 #### Gaussian approximation of proposal using Taylor Expansion up to 2.
 
 propvarfun2 <- function(x, mu, od, data, theta, tausq){
   atmp = ( ((x-xi2(data,mu,od))/sqrt(tausq))*dnorm((xi2(data,mu,od)-x)/sqrt(tausq)) - 
              ((x-xi2(data-1,mu,od))/sqrt(tausq))*dnorm((xi2(data-1,mu,od)-x)/sqrt(tausq)) )/
     (pnorm((xi2(data,mu,od)-x)/sqrt(tausq)) - pnorm((xi2(data-1,mu,od)-x)/sqrt(tausq)))
   
   btmp = (( dnorm((xi2(data-1,mu,od)-x)/sqrt(tausq))-dnorm((xi2(data,mu,od)-x)/sqrt(tausq)) )/
             ( pnorm((xi2(data,mu,od)-x)/sqrt(tausq))-pnorm((xi2(data-1,mu,od)-x)/sqrt(tausq))))^2
   
   c2 = (atmp-btmp)/tausq
   corr1 = solve((1-tausq)*exp(-distwood/theta))
   cov2 = (solve(corr1-diag(c2)))
   return(cov2)
 }
 
 #### Initial Values
 beta0 = estzip00$estimate[1:2]
 od0 = estzip00$estimate[3]
 theta0 = estzip00$estimate[4]
 tausq0 = estzip00$estimate[5]
 
 #### Initial Value of gamma0
 z01 = qnorm(pzip(mapley-1, od = od0, mu = exp(beta0[1]+beta0[2]*hickx)))
 z02 = qnorm(pzip(mapley, od = od0, mu = exp(beta0[1]+beta0[2]*hickx)))
 z0 = (dnorm(z01)-dnorm(z02))/(pnorm(z02)-pnorm(z01))
 
 #### Posterior mode
 tmpZip0 = spg(par = z0, fn = llik11, mu = exp(beta0[1]+beta0[2]*hickx), od = od0, data = mapley,
               theta = theta0, tausq = tausq0, 
               control = c(maxit = 10000, ftol = 1e-5, gtol = 1e-5,trace = F))
 
 #### Posterior Var if using normal approximation with Taylor Expansion up to 2
 propmuzip = tmpZip0$par
 propvarzip = propvarfun2(x = tmpZip0$par, mu = exp(beta0[1]+beta0[2]*hickx), od = od0, data = mapley,
                          theta = theta0, tausq = tausq0)
 
 
 #####################################################################################################################
 ############# SET UP THE CHAIN
 #####################################################################################################################
 nit <- 5000 #nit <- 15000 
 k = 40
 theta <- rep(0, nit); tausq <- rep(0, nit); od <- rep(0, nit)
 beta <- matrix(0, nrow = nit, ncol = 2)
 gamma <- array(0, dim = c(256, k, nit))
 mu <- matrix(0, nrow = nit, ncol = 256)
 
 gammazip0 <- t(mvrnorm(k, mu = propmuzip, Sigma = propvarzip))
 
 # Initialize of all parameters including the latent fields with k copies
 theta[1] <- theta0; tausq[1] <- tausq0; od[1] <- od0
 beta[1,] <- beta0;  gamma[,,1] <- gammazip0
 mu[1,] <- exp(beta0[1] + beta0[2]*hickx)
 
 ####
 #### Start the Chain
 #### Two Blocks Langevin Updates
 
 counter1 <- 0; counter2 <- 0
 
 #########################################################################################################################
 ######################################  The Proposals when K = 40 Are As Following  #####################################
 #########################################################################################################################
 
 # Langevin stepsize of beta1, beta2, gamma
 h1 = 0.000085; h2 = 0.0000018; h3 = 0.088 # maybe try 0.09
 # Langevin stepsize of beta1, beta2, gamma
 dd = 0.0007; tt = 0.00066; gg = 0.0011  
 
 # Langevin stepsize of logit(tausq)
 #trans.nod <- seq(500,nit,by = 500)  
 # Where Adaptive may Take Place
 
 for (i in 1:(nit-1)){
   
   # For RWP we can use this 
   #  if (i %in% trans.nod){
   #    step.beta = as.vector(cov(beta[(i-499):i,]))*(2.38/sqrt(2))^2
   #    step.theta = sd(log(theta[(i-499):i]))*2.38/sqrt(2)
   #    step.other.sd[2] = sd(logit(tausq[(i-499):i]))*2.38/sqrt(2)
   #  }
   #### Calculate Gradient of beta1 and beta2 and
   #### Set Langevin Updates
   
   grad11 = grad(function(x) {betapost12(x, od = od[i], data = mapley, tausq = tausq[i], 
                                         gammas = gamma[,,i], beta2 = beta[i,2])},beta[i,1]) - (beta[i,1])/30
   
   bmunew1 = beta[i,1] + (h1/2)*grad11
   betanew1 <- rnorm(1, mean = bmunew1, sd = sqrt(h1))
   grad12 = grad(function(x) {betapost12(x, od = od[i], data = mapley, tausq = tausq[i], 
                                         gammas = gamma[,,i], beta2 = beta[i,2])},betanew1) - (betanew1)/30
   bmu1back <- betanew1+(h1/2)*grad12
   grad21 =  grad(function(x) {betapost22(x, od = od[i], data = mapley, tausq = tausq[i], 
                                          gammas = gamma[,,i], beta1 = beta[i,1])},beta[i,2]) - (beta[i,2])/30
   
   bmunew2 = beta[i,2] + (h2/2)*grad21
   betanew2 <- rnorm(1, mean = bmunew2, sd = sqrt(h2))
   munew <- exp(betanew1+betanew2*hickx)
   
   grad22 =  grad(function(x) {betapost22(x, od = od[i], data = mapley, tausq = tausq[i], 
                                          gammas = gamma[,,i], beta1 = beta[i,1])},betanew2) - (betanew2+0.1)/30
   bmu2back <- betanew2 + (h2/2)*grad22
   
   gradlod <- grad(function(x) {odpost(x, beta1 = beta[i,1], beta2 = beta[i,2], 
                                       data = mapley, tausq = tausq[i],  gammas = gamma[,,i])},od[i])*od[i]
   lod <- log(od[i])
   mulodnew <- lod+0.5*dd*gradlod
   lodnew <- rnorm(1, mean = mulodnew, sd = sqrt(dd))
   odnew <- exp(lodnew)
   
   gradlodnew <- grad(function(x) {odpost(x, beta1 = beta[i,1], beta2 = beta[i,2], 
                                          data = mapley, tausq = tausq[i],  gammas = gamma[,,i])},odnew)*odnew
   mulodback <- lodnew+0.5*dd*gradlodnew
   
   logratio1 <- min(0, (dnorm(betanew1, mean = 0, sd = sqrt(30), log = T)+
                          dnorm(betanew2, mean = 0, sd = sqrt(30), log = T)-
                          dnorm(beta[i,1], mean = 0, sd = sqrt(30), log = T)-
                          dnorm(beta[i,2], mean = 0, sd = sqrt(30), log = T)+
                          dlnorm(odnew, meanlog = 0, sdlog = 3, log = T) - 
                          dlnorm(od[i], meanlog = 0, sdlog = 3, log = T)+
                          pzip1(mapley, tausq[i], gamma[,,i], munew, odnew)- 
                          pzip1(mapley, tausq[i], gamma[,,i], mu[i, ], od[i])+
                          dnorm(beta[i, 1], mean = bmunew1, sd = sqrt(h1), log = T) -
                          dnorm(betanew1, mean = bmu1back, sd = sqrt(h1), log = T) +
                          dnorm(beta[i, 2], mean = bmunew2, sd = sqrt(h2), log = T) -
                          dnorm(betanew2, mean = bmu2back, sd = sqrt(h2), log = T) +
                          dnorm(lod, mean = mulodnew, sd = sqrt(dd), log = T) -
                          dnorm(lodnew, mean = mulodback, sd = sqrt(dd), log = T)+
                          log(odnew)-log(od[i])
   ))
   
   if (runif(1) <= exp(logratio1)){
     beta[i+1,1] = betanew1; beta[i+1,2] = betanew2; od[i+1] = odnew; mu[i+1, ] = munew; counter1 = counter1+1
   } else {
     beta[i+1,1] = beta[i,1]; beta[i+1,2] = beta[i,2]; od[i+1] = od[i]; 
     mu[i+1,] = mu[i,]; counter1 = counter1 
   }
   
   #### Second Block: Update latent fields and the dependence parameters  
   #### For tausq: using Langevin Update to logit(tausq)
   #### For theta: using rw to log(theta)
   #### Tentatively work fine, may be improved anyway.
   
   gradltheta <- grad(function(x) {post.theta(x, data = mapley, gammas = gamma[,,i], 
                                              tausq = tausq[i])},theta[i])*theta[i]
   ltheta <- log(theta[i])
   multhetanew <- ltheta+0.5*tt*gradltheta
   lthetanew <- rnorm(1, mean = multhetanew, sd = sqrt(tt))
   thetanew <- exp(lthetanew)
   
   gradlthetanew <- grad(function(x) {post.theta(x, data = mapley, gammas = gamma[,,i], 
                                                 tausq = tausq[i])},thetanew)*thetanew
   multhetaback <- lthetanew+0.5*tt*gradlthetanew
  
   phi <- logit(tausq[i])
   muphinew <- phi+0.5*gg*grad.logit.tausq2(mapley, od[i], gamma[,,i], mu[i+1, ], theta[i], tausq[i])
   phinew <- rnorm(1, mean = muphinew, sd = sqrt(gg))
   
   tausqnew <- inv.logit(phinew)
   tausqnew <- ifelse(tausqnew > 0.99, 0.99, tausqnew) 
   tausqnew <- ifelse(tausqnew < 0.0001, 0.0001, tausqnew) 
   
   muphiback <- phinew + 0.5*gg*grad.logit.tausq2(mapley, od[i], gamma[,,i], mu[i+1, ], theta[i], tausqnew)
   grad.mat1 <- grad.joint2(x = gamma[,,i], mu = mu[i+1,], od = od[i+1], 
                            data = mapley, theta = theta[i], tausq = tausq[i])
   mu.gammanew = gamma[,,i] + (h3/2)*propvarzip%*%grad.mat1
   gammanew = matmvnrmCpp(mu = mu.gammanew, Sigma = h3*propvarzip)
   grad.mat2 <- grad.joint2(x = gammanew, mu = mu[i+1,], od = od[i+1], 
                            data = mapley, theta = theta[i], tausq = tausq[i])
   mu.gammaback = gammanew + (h3/2)*propvarzip%*%grad.mat2
   
   logratio2 <- min(0, (    pzip1(mapley,tausqnew, gammanew, mu[i+1,], od[i+1])- 
                              pzip1(mapley, tausq[i], gamma[,,i], mu[i+1, ], od[i+1])+
                              dlmvnormx(t(gammanew),mean=rep(0,256),sigma = (1-tausqnew)*exp(-distwood/thetanew)) -
                              dlmvnormx(t(gamma[,,i]),mean=rep(0,256),sigma = (1-tausq[i])*exp(-distwood/theta[i])) +
                              dnorm(phi, mean = muphinew, sd = sqrt(gg), log = T) -
                              dnorm(phinew, mean = muphiback, sd = sqrt(gg), log = T) +
                              dnorm(ltheta, mean = multhetanew, sd = sqrt(tt), log = T) -
                              dnorm(lthetanew, mean = multhetaback, sd = sqrt(tt), log = T)+
                              log(tausqnew)+log(1-tausqnew)-log(tausq[i])-log(1-tausq[i])+log(thetanew)-log(theta[i]) + 
                              #The Jacobian
                              dlnorm(theta[i], meanlog = -2, sdlog = 1, log = T)-dlnorm(thetanew, meanlog = -2, sdlog= 1, log = T) + 
                              #Arbitrary Prior
                              dlmvnmatCpp(gamma[,,i], mu = mu.gammaback, Sigma = h3*propvarzip)-
                              dlmvnmatCpp(gammanew, mu = mu.gammanew, Sigma = h3*propvarzip)
   ))
   if (runif(1) <= exp(logratio2)){
     tausq[i+1] = tausqnew; theta[i+1] = thetanew; gamma[,,i+1] = gammanew; counter2 = counter2+1
   } else {
     theta[i+1] = theta[i]; tausq[i+1] = tausq[i]; gamma[,,i+1] = gamma[,,i]; counter2 = counter2
   }
   
   #### Save too many gamma will slow down your computer and eat lots of memory  !!!
   #### Delete one by one will resolve the issue and increase speed.
   #####
 }
 
 # > mean(beta[5000:nit,1]) 
 # [1] 0.8950802
 # > mean(beta[5000:nit,2]) 
 # [1] -0.1553258
 # > mean(od[5000:nit]) 
 # [1] 0.3554797
 # > mean(theta[5000:nit]) 
 # [1] 0.08248159
 # > mean(tausq[5000:nit]) 
 # [1] 0.4084383
 #> sqrt(k*(var(beta[5000:nit,1]))) 
 #[1] 0.1312127
 #> sqrt(k*var(beta[5000:nit,2]))  
 #1] 0.02629007
 #> sqrt(k*var(od[5000:nit]))   
 #[1] 0.1077338
 #> sqrt(k*var(theta[5000:nit]))  
 #[1] 0.02064321
 #> sqrt(k*var(tausq[5000:nit]))
 #[1] 0.1155504
 
 
##########################################################################################################
 ################# Using GB (Nikoloulopoulos 2013) method 
 #########################################################################################################
  directlik2 <- function(v)
 {
   #v = c(1.1,-0.15,0.08,0.5)
   beta0 = v[1]; beta1 = v[2]; od = v[3]; theta = v[4]; tausq = v[5]
   mu = exp(beta0 + beta1*hickx)
   
   lower = qnorm(pzip(q = mapley-1, od = od, mu = mu))
   upper = qnorm(pzip(q = mapley, od = od, mu = mu))
   
   R = (1-tausq)*exp(-distwood/theta)+(distwood==0)*tausq
   
   set.seed(1234)
   lik = -log(pmvnorm(lower = lower, upper = upper, mean = rep(0, length(lower)), sigma = R)[1])
   
   lik = ifelse(lik == Inf, -log(.Machine$double.xmin), lik)
   lik = ifelse(lik == -Inf, log(.Machine$double.xmin), lik)
   return(lik)
 }
 
 gbresult <- try(optim(par = c(1.2,-0.2,0.07,0.5), fn = directlik2, gr = NULL, method = c("L-BFGS-B"),
                       lower = c(-3, -3, 0.01,0.01, 0.01), upper = c(5, 5, 4, 5, 0.95), hessian = T), silent = T)
 
 
 gbresult.gbzip <- optim(par = c(0.9, -0.15, 0.35, 0.08, 0.41), fn = directlik2, method = "L-BFGS-B", 
                         lower = c(-Inf, -Inf, 0, 0, 0), upper = c(Inf, Inf, Inf, Inf, 0.9), hessian = T)

 
 
 
 
 ################################################################################################################
 ################################################################################################################
 ############################ Analyze the Lansing Woods Data with NB Marginal ##################################
 ################################################################################################################
 ################################################################################################################
 ### GHK estimates via Package gcKrig
 estnbGHK <- mlegc(y = LansingTrees$maple, x = LansingTrees$hickory, 
                    locs = cbind(LansingTrees$Easting, LansingTrees$Northing), 
                    marginal = negbin.gc(link = 'log'), 
                    corr = matern.gc(nugget = TRUE))
 
 ### GHK estimaes with gcmr; closed to the gcKrig
 estnbgcmr <- gcmr(mapley~hickx, marginal = negbin.marg(link='log'), cormat = matern_nugget.cormat(distwood),
                   options = gcmr.options(seed = 1234, nrep = c(200,5000)))
 
 ### Compute Residuals
 estnbgcmr$estimate <- estnbGHK$MLE
 residNB <- gcmr::residuals(estnbgcmr, type = "conditional", method = "random")
 qqnorm(residNB, main = NULL)
 qqline(residNB, col = 2, lwd = 2)
 #dev.print(file = "ResidualQQ.eps")
 plot(residNB, xlab = "Location Index", ylab = "Residual", main = NULL)
 abline(h = 0, col = 2, lwd = 1.8)
 dev.print(file = "Residual.eps")
 hist(residNB, density=20, breaks=20, prob=TRUE, 
      xlab="Residual", ylim=c(0, 0.5), main = NULL, xlim = c(-3,3))
 curve(dnorm(x, mean=0, sd=1), 
       col="darkblue", lwd=2, add=TRUE, yaxt="n")
 #dev.print(file = "ResidualHist.eps")
 acf(residNB, main = NA)
 #dev.print(file = "ResidualACF.eps")

 ### Run 500 with Randomized quantile residuals; test with Shapiro Wilk
 set.seed(1234)
 seedresid <- sample(1:10000, 500, replace = FALSE) 
 residualp <- rep(0,500)
 for(i in 1:500){
 residNB <- residuals(estnbgcmr, type = "conditional", method = "random")
 residualp[i] <- shapiro.test(residNB)$p.value
 }
 hist(residualp, density=20, breaks=20, prob=TRUE, 
      xlab="p-value", ylim=c(0, 4), main = NULL)
 #dev.print(file = "pHist.eps")
 
 #####################################################################################################################
 ##################### FOR SIMULATED DATA ############################
 #####################################################################################################################
 
 xloc=rep(seq(0,1,by=0.1),11)                  
 yloc=rep(seq(0,1,by=0.1),each=11)             
 mat=matrix(0,121,121)                        
 mat=as.matrix(dist(cbind(xloc,yloc),method='euclidean'))
 
 set.seed(529529)
 seedsim <- sample(99999999,1000)
 
 data.nbsim1 <- matrix(0, ncol = 121, nrow = 1000)
 for (i in 1:1000){
   set.seed(seedsim[i])
   data.nbsim[i,] <- qnbinom(pnorm(mvrnorm(1,mu=rep(0,121),Sigma = (1-0.25)*exp(-mat/0.1)+(mat==0)*0.25)),
                             size= rep(0.5, 121), mu = exp(0.5+ 0.5*xloc+ yloc))
 }
 
 dlmvnormx <- function(x,mean,sigma){
   llik <- dmvnorm(x, mean = mean, sigma = sigma, log = T)
   llik <-  ifelse(llik==-Inf, -.Machine$double.xmax, llik)
   llik <- ifelse(llik== Inf, .Machine$double.xmax, llik)
   llik[is.nan(llik)] == log(.Machine$double.xmin)
   return(sum(llik))
 }
 
 #### Log likelihood of: prod(hi(yi|latent, parameters)) with cloned latents
 pnb1 <- function(data, tausq, gammas, mu, od){
   prob = (pnorm((qnorm(pnbinom(q = data, size = 1/od, mu = mu))-gammas)/sqrt(tausq))-
             pnorm((qnorm(pnbinom(q = data-1, size = 1/od, mu = mu))-gammas)/sqrt(tausq)))
   prob2 = ifelse(prob==0, .Machine$double.eps, prob)
   sumloglike2 = sum(log(prob2))
   return(sumloglike2)
 }
 
 
 xi <- function(x, mu, od) {
   out = qnorm(pnbinom(x, size = 1/od, mu = mu))
   out[is.infinite(out)] = 0
   return(out)
 }
 
 grad.joint <- function(x, mu, od, data, theta, tausq){
   num = dnorm((qnorm(pnbinom(data-1, size = 1/od, mu = mu))-x)/sqrt(tausq)) - 
     dnorm((qnorm(pnbinom(data, size = 1/od, mu = mu))-x)/sqrt(tausq))
   denum = pnorm((qnorm(pnbinom(data, size = 1/od, mu = mu))-x)/sqrt(tausq)) - 
     pnorm((qnorm(pnbinom(data-1, size = 1/od, mu = mu))-x)/sqrt(tausq))
   denum[denum ==0] = .Machine$double.eps
   gradmvn = solve( (1-tausq)*exp(-mat/theta))%*%(-x) 
   return(gradmvn + num/(denum*sqrt(tausq)))
 }
 
 est0 <- mlegc(y = data.nbsim[1,], x = cbind(xloc, yloc), locs = cbind(xloc, yloc), 
               marginal = negbin.gc(link = 'log'), corr = matern.gc(nugget = TRUE))
 
 est0 <- est0$MLE
 nit <- 12000; k = 30
 theta <- rep(0, nit); beta <- matrix(0, nrow = nit, ncol = 3)
 od <- rep(0,nit); tausq <- rep(0, nit)
 gamma <- array(0, dim = c(121, k, nit))
 mu <- matrix(0, nrow = nit, ncol = 121)
 
 beta[1,] = est0[1:3]
 od[1] = est0[4]
 theta[1] = est0[5]
 tausq[1] = ifelse(est0[6] == 1, 0.8, est0[6])
 mu[1,] = exp(beta[1,1] + beta[1,2]*xloc + beta[1,3]*yloc)
 gamma[,,1] = t(mvrnorm(k, mu = rep(0,121), Sigma = (1-tausq[1])*exp(-mat/theta[1])))
 
 #### Compute log-likelihood of mvn matrix, with different mean vector. 
 cppFunction(depends = 'RcppArmadillo',   code = '
             double dlmvnmatCpp(arma::mat x, arma::mat mu, arma::mat Sigma ) { 
             const double log2pi = std::log(2.0 * M_PI);
             int n = x.n_cols; int xdim = x.n_rows;
             arma::vec out(n);
             double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
             arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
             double rootisum = arma::sum(log(rooti.diag()));
             for (int i=0; i < n; i++) {
             arma::vec z = rooti * (x.col(i) - mu.col(i));
             out(i) = constants - 0.5 * arma::sum(z%z) + rootisum;
             if(out(i) == 0 ){
             out(i) =  std::numeric_limits<float>::min() ;
             }
             }
             return arma::sum(out);
             }')

 cppFunction(depends = 'RcppArmadillo',   code = '
             arma::mat matmvnrmCpp (arma::mat mu, arma::mat Sigma){
             int k = mu.n_cols;
             int n = mu.n_rows;
             arma::mat Y = arma::randn(k, n);
             return (mu.t() + Y*arma::chol(Sigma)).t(); 
             }')

 counter1 <- 0; counter2 <- 0
 hh = c(0.0024, 0.0048, 0.0054)
 dd = 0.003
 tt = 0.00027
 gg = 0.0016
 ll = 0.03 
 
 for (i in 1:(nit-1)){
   
   betanew <- mvrnorm(mu = beta[i,], Sigma = diag(hh))
   munew <- exp(betanew[1] + betanew[2]*xloc + betanew[3]*yloc)
   lod <- log(od[i])
   lodnew <- rnorm(1, mean = lod, sd = sqrt(dd))
   odnew <- exp(lodnew)
   
   logratio1 <- min(0, ( dmvnorm(betanew, mean = c(0.5,0.5,1), sigma = diag(30, 3), log = T) - 
                           dmvnorm(beta[i,], mean = c(0.5,0.5,1), sigma = diag(30, 3), log = T)+
                           dlnorm(odnew, meanlog = 0, sdlog = 3, log = T) - 
                           dlnorm(od[i], meanlog = 0, sdlog = 3, log = T)+
                           pnb1(data.nbsim[1,], tausq[i], gamma[,,i], munew, odnew)- 
                           pnb1(data.nbsim[1,], tausq[i], gamma[,,i], mu[i, ], od[i])+
                           log(odnew)-log(od[i])
   ))
   
   if (runif(1) <= exp(logratio1)){
     beta[i+1, ] = betanew; od[i+1] = odnew; mu[i+1, ] = munew; counter1 = counter1+1
   } else {
     beta[i+1, ] = beta[i, ]; od[i+1] = od[i]; mu[i+1,] = mu[i,]; counter1 = counter1 
   }
   
   #### Second Block: Update latent fields and the dependence parameters  
   ltheta <- log(theta[i])
   lthetanew <- rnorm(1, mean = ltheta, sd = sqrt(tt))
   thetanew <- exp(lthetanew)  
   
   phi <- log(tausq[i]/(1-tausq[i]))
   phinew <- rnorm(1, mean = phi, sd = sqrt(gg))
   tausqnew <- exp(phinew)/(1+exp(phinew))
   
   #### UPDATE LATENT FIELD      
   
   grad.mat1 <- grad.joint(x = gamma[,,i], mu = mu[i+1,], od = od[i+1], 
                           data = data.nbsim[1,], theta = theta[i], tausq = tausq[i])
   mu.gammanew = gamma[,,i] + (ll/2)*grad.mat1
   
   gammanew = matmvnrmCpp(mu = mu.gammanew, Sigma = diag(ll,121))
   grad.mat2 <- grad.joint(x = gammanew, mu = mu[i+1,], od = od[i+1], 
                           data = data.nbsim[1,], theta = theta[i], tausq = tausq[i])
   mu.gammaback = gammanew + (ll/2)*grad.mat2
   
   logratio2 <- min(0, (pnb1(data.nbsim[1,],tausqnew, gammanew, mu[i+1,], od[i+1])- 
                          pnb1(data.nbsim[1,], tausq[i], gamma[,,i], mu[i+1, ], od[i+1])+
                          dlmvnormx(t(gammanew),mean=rep(0,121),sigma = (1-tausqnew)*exp(-mat/thetanew)) -
                          dlmvnormx(t(gamma[,,i]),mean=rep(0,121),sigma = (1-tausq[i])*exp(-mat/theta[i])) +
                          log(tausqnew)+log(1-tausqnew)-log(tausq[i])-log(1-tausq[i])+log(thetanew)-log(theta[i]) +
                          dlnorm(theta[i], meanlog = -2, sdlog = 1, log = T) - 
                          dlnorm(thetanew, meanlog = -2, sdlog= 1, log = T) + 
                          dlmvnmatCpp(gamma[,,i], mu = mu.gammaback, Sigma = diag(ll,121))-
                          dlmvnmatCpp(gammanew, mu = mu.gammanew, Sigma = diag(ll,121))
   ))
   
   if (runif(1) <= exp(logratio2)){
     tausq[i+1] = tausqnew; theta[i+1] = thetanew; gamma[,,i+1] = gammanew; counter2 = counter2+1
   } else {
     theta[i+1] = theta[i]; tausq[i+1] = tausq[i]
     gamma[,,i+1] = gamma[,,i]; counter2 = counter2
   }
 }
 