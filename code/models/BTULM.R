fit_dep_gibbs <- function(n_gibbs=2000, n_burn=500, sigma2_beta=10000,
                          Psi,X, Y, weights, weeks, areas){
start.time <- Sys.time()
D <- ncol(Psi)
n_weeks <- max(weeks) - min(weeks) + 1
p <- ncol(X)
N <- nrow(X) #number of responses (including repeat respondents)

I_p <- Diagonal(p)
I_D <- Diagonal(D)
zero_p <- rep(0, p)
zero_D <- rep(0, D)

###
### Setup storage
###
n_keep <- n_gibbs - n_burn
eta.post <- array(NA, dim = c(n_keep, D, n_weeks))
betas.post <- matrix(NA, n_keep, p)
omegas.post <- array(NA, dim=c(n_keep, N))
phi.post <- sigma2_eta.post <- sigma2_eta_1.post <-rep(NA, n_keep)

###
### Set initial values
###
sigma2_eta <- 1
sigma_eta <- sqrt(sigma2_eta)
sigma2_eta_1 <- 1
phi <- .5
beta <- glm(formula=as.formula(paste("RESPONSE" ,"~", pop_covars)),
             family=binomial,
             data=HPS_pop_long) %>%
         coefficients %>% 
         as.vector

eta <- Matrix(0, nrow=D, ncol=n_weeks)
for(t in 1:n_weeks){
  eta[,t] <- rmvn(1, rep(0, D), diag(rep(.1, D)))
}

#IG parameters
a <- .1
b <- .1

###
### Start sampling
###

kappa <- weights*(Y - .5)

pb <- txtProgressBar(min=0, max=n_gibbs, style=3)
for(iter in 1:n_gibbs){
  Psi_eta <- (Psi %*% eta)[cbind(seq_along(weeks), weeks)]
  #sample latent PG variable
  omega <- rpg.gamma(num=N, h=weights, 
               z=as.numeric(X%*%beta + Psi_eta))

  ##sample sigma_eta_1
  sigma2_eta_1 <-  rinvgamma(n=1, shape=a + D/2,
                              scale=b + t(eta[,1])%*%eta[,1]/2)

  #sample sigma^2_eta
  sigma2_eta_shape <- a + D*(n_weeks-1)/2
  eta_diff <- eta[,2:n_weeks] - phi * eta[,1:(n_weeks-1)]
  sigma2_eta_scale <- b + sum(eta_diff * eta_diff)/2

  sigma2_eta <- rinvgamma(n=1, shape=sigma2_eta_shape,
                               scale=sigma2_eta_scale)
  sigma_eta <- sqrt(sigma2_eta)
 
  #sample phi
  denom <- sum(eta[,1:(n_weeks-1)] * eta[,1:(n_weeks-1)]) 
  mu_phi <- sum(eta[,2:n_weeks] * eta[,1:(n_weeks-1)])/denom
  sigma_phi <- sigma_eta/sqrt(denom)
  phi <- rtrunc(n=1, spec="norm", a=-1, b=1, mean=mu_phi, sd=sigma_phi)

  #sample beta
  Omega <- Diagonal(x=omega)
  time1<-Sys.time()
  sqrt_XOmega <- sqrt(Omega)%*%X
  Prec_beta <- t(sqrt_XOmega)%*%(sqrt_XOmega) + I_p/sigma2_beta
  mu_beta <- t(X)%*%Omega%*%(kappa/omega-Psi_eta)

  U_beta <- base::chol(Prec_beta)
  tmp_norm <- rnorm(nrow(Prec_beta))
  beta <- backsolve(U_beta, backsolve(U_beta, mu_beta, transpose=TRUE) + tmp_norm)

  #sample eta[,1]
  ids_1 <- which(weeks==1)
  Psi_1 <- Psi[ids_1,]
  omega_1 <- omega[ids_1]
  Omega_1 <- Diagonal(x=omega_1)
  X_1 <- X[ids_1,]
  kappa_1 <- kappa[ids_1]

  Prec_eta_1 <- t(Psi_1)%*%Omega_1%*%Psi_1 + 
                         (1/sigma2_eta_1 + phi^2/sigma2_eta)*I_D

  mu_eta_1 <- t(Psi_1)%*%Omega_1%*%(kappa_1/omega_1 - X_1%*%beta) + phi/sigma2_eta*eta[,2]

  U_eta_1 <- base::chol(Prec_eta_1)
  tmp_norm <- rnorm(nrow(Prec_eta_1))
  eta[,1] <- backsolve(U_eta_1, 
                       backsolve(U_eta_1, mu_eta_1, 
                                 transpose=TRUE) + tmp_norm)

  ##sample eta[,2:(n_weeks-1)]
  for(t in 2:(n_weeks-1)){
    ids_t <- which(weeks==t)
    Psi_t <- Psi[ids_t,]
    omega_t <- omega[ids_t]
    Omega_t <- Diagonal(x=omega_t)
    X_t <- X[ids_t,]
    kappa_t <- kappa[ids_t]

    Prec_eta_t <- t(Psi_t)%*%Omega_t%*%Psi_t + (1+phi^2)/sigma2_eta * I_D
    mu_eta_t <- t(Psi_t)%*%Omega_t%*%(kappa_t/omega_t - X_t%*%beta) + phi/sigma2_eta*(eta[,(t-1)] + eta[,(t+1)])

    U_eta_t <- base::chol(Prec_eta_t)
    tmp_norm <- rnorm(nrow(Prec_eta_t))
    eta[,t] <- backsolve(U_eta_t, backsolve(U_eta_t, mu_eta_t, transpose=TRUE) + tmp_norm)
  }

  #sample eta[,n_weeks]
  ids_T <- which(weeks==n_weeks)
  Psi_T <- Psi[ids_T,]
  omega_T <- omega[ids_T]
  Omega_T <- Diagonal(x=omega_T)
  X_T <- X[ids_T,]
  kappa_T <- kappa[ids_T]

  Prec_eta_T <- t(Psi_T)%*%Omega_T%*%Psi_T + I_D/sigma2_eta
  mu_eta_T <- t(Psi_T)%*%Omega_T%*%(kappa_T/omega_T - X_T%*%beta) + phi/sigma2_eta*eta[,(n_weeks-1)]
  
  U_eta_T <- base::chol(Prec_eta_T)
  tmp_norm <- rnorm(nrow(Prec_eta_T))
  eta[,n_weeks] <- backsolve(U_eta_T, backsolve(U_eta_T, mu_eta_T, transpose=TRUE) + tmp_norm)
  
  if(iter > n_burn){#discard burn-in steps
   index <- iter-n_burn
   eta.post[index,,] <- as.matrix(eta)
   betas.post[index,] <- as.vector(beta)
   phi.post[index] <- phi
   sigma2_eta.post[index] <- sigma2_eta
   sigma2_eta_1.post[index] <- sigma2_eta_1
   }
 setTxtProgressBar(pb, iter)
 }
 colnames(betas.post) <- colnames(X)
 runtime <- as.numeric(difftime(Sys.time(), start.time), units="secs")
 posteriors <- list(eta=eta.post, betas=betas.post, phi=phi.post,
                    sigma2_eta=sigma2_eta.post, #omegas=omegas.post,
                    sigma2_eta_1=sigma2_eta_1.post, runtime=runtime)
 return(posteriors)
}
