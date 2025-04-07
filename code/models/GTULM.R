fit_dep_gibbs <- function(n_gibbs=2000, n_burn=500,
                      y_1, y_2, N_1, N_2, 
                      n_weeks, n_areas, X_1, X_2,  
                      areas_1, areas_2, weeks_1, 
                      weeks_2, weights_1, weights_2){

  weeks_1 <- weeks_1 - min(weeks_1) + 1
  weeks_2 <- weeks_2 - min(weeks_2) + 2

  X_1 <- Matrix(X_1)
  X_2 <- Matrix(X_2)
  ###
  ### SET USEFUL VARIABLES
  ###
  N <- sum(N_1) + sum(N_2)
  D <- n_areas 
  #T <- n_weeks
  p <- ncol(X_1)
  W_1 <- Diagonal(x=weights_1)
  W_2 <- Diagonal(x=weights_2)
  
  #create incidence matrices
  Psi_1 <- Matrix(model.matrix(~factor(areas_1)-1))
  Psi_2 <- Matrix(model.matrix(~factor(areas_2)-1))
  
  sigma2_beta <- 10000 
  I_p <- Diagonal(p)
  I_D <- Diagonal(D)
  zero_p <- rep(0, p)
  zero_D <- rep(0, D)

  ###
  ### SETUP STORAGE
  ### 
  n_keep <- n_gibbs - n_burn
  eta.post <- array(NA, dim = c(n_keep, n_areas, n_weeks))
  betas.post <- matrix(NA, n_keep, p)
  rho.post <- phi.post <- sigma2.post <- sigma2_eta.post <- sigma2_eta_1.post <-rep(NA, n_keep)
  
  ###
  ### SET INITIAL VALUES
  ###
  sigma2_eta <- 1
  sigma_eta <- sqrt(sigma2_eta)
  sigma2_eta_1 <- 1
  sigma_eta_1 <- sqrt(sigma2_eta_1)
  phi <- .5#runif(n=1, -1, 1)
  beta <- t(rmvn(n=1, zero_p, Sigma=diag(1, nrow=p)))
  rho <- .5 
  rho.num_accept <- 0
  
  eta <- Matrix(0, nrow=D, ncol=n_weeks)
  eta[, 1] <- rmvn(n=1, zero_D, Sigma=diag(sigma2_eta_1, nrow=D))
  for(t in 2:n_weeks){
    eta[, t] <- rmvn(n=1, phi*eta[, (t-1)], Sigma=diag(sigma2_eta, nrow=D))
  }
  
  ###
  ### START SAMPLING
  ###
  for(iter in 1:n_gibbs){
    if(iter %% 50 == 0){print(paste0("Gibbs iteration #",iter))}
    ##sample phi
    denom <- sum(eta[,1:(n_weeks-1)] * eta[,1:(n_weeks-1)]) #quicker way to get trace
    mu_phi <- sum(eta[,2:n_weeks] * eta[,1:(n_weeks-1)])/denom
    sigma_phi <- sigma_eta/sqrt(denom)
    phi <- rtrunc(n=1, spec="norm", a=-1, b=1, mean=mu_phi, sd=sigma_phi)
    
    ##sample sigma_eta
    sigma2_eta_shape <- .1 + .5 * D * (n_weeks-1)
    eta_diff <- eta[,2:n_weeks] - phi * eta[,1:(n_weeks-1)]
    sigma2_eta_scale <- .1 + .5 * sum(eta_diff * eta_diff)  #quicker way to get trace
    sigma2_eta <- rinvgamma(n=1, shape=sigma2_eta_shape, 
                                 scale=sigma2_eta_scale)
    sigma_eta <- sqrt(sigma2_eta)

    ##sample sigma_eta_1
    sigma2_eta_1 <- rinvgamma(n=1, shape=.1 + D/2,
                                 scale=.1 + .5*t(eta[,1])%*%eta[,1])
    sigma_eta_1 <- sqrt(sigma2_eta_1)
  
    #sample sigma
    #using this trick for faster indexing: https://stackoverflow.com/a/25584122
    # we extract the vector of etas corresponding to the random effects for 
    # all first time respondents, then for all repeat respondents (both initial
    # and follow-up response)
    eta_f <- (Psi_1 %*% eta)[cbind(seq_along(weeks_1), weeks_1)]
    eta_r2 <- (Psi_2 %*% eta)[cbind(seq_along(weeks_2), weeks_2)]
    eta_r1 <- (Psi_2 %*% eta)[cbind(seq_along(weeks_2), (weeks_2-1))]

    y_tilde_f <- y_1 - X_1%*%beta - eta_f
    y_tilde_r2 <- y_2[,2] - X_2%*%beta - eta_r2
    y_tilde_r1 <- y_2[,1] - X_2%*%beta - eta_r1
    y_tilde_r_diff <- y_tilde_r2-rho*y_tilde_r1
    sigma2_shape <- 1 + N/2
    sigma2_scale <- 1 + 1/2 *( (t(y_tilde_f) %*% W_1 %*% y_tilde_f) + 
  		  1/(1-rho^2) * (t(y_tilde_r_diff) %*% W_2 %*% y_tilde_r_diff))
  
    sigma2 <- rinvgamma(n=1, shape=sigma2_shape, 
                             scale=as.numeric(sigma2_scale))
    sigma <- sqrt(sigma2)
    
    #sample rho with Metropolis step
    spread <- min(abs(-1-rho), abs(1-rho), 0.05)
    rho_proposal <- runif(n=1, rho-spread, rho+spread)
    
    rn <- sum(weights_2*dnorm(y_2[,2], 
                              mean=as.vector((X_2%*%beta + eta_r2) + rho_proposal*(y_tilde_r1)),
                              sd=sqrt(sigma2*(1-rho_proposal^2)),
                              log=T))
    rd <- sum(weights_2*dnorm(y_2[,2], 
                              mean=as.vector((X_2%*%beta + eta_r2) + rho*(y_tilde_r1)),
                              sd=sqrt(sigma2*(1-rho^2)),
                              log=T))
    ratio <- exp(rn-rd)
    threshold <- runif(1, 0, 1)
    if((ratio) > threshold){
      rho <- rho_proposal
      if(iter>n_burn){
        rho.num_accept <- rho.num_accept + 1
      }
    }
    
    #sample beta
    Xft_Wf_Xf <- t(X_1) %*% W_1 %*% X_1 
    Xrt_Wr_Xr <- t(X_2) %*% W_2 %*% X_2 
 
    Prec_beta <- 1/sigma2 * (Xft_Wf_Xf - (1-rho)/(1+rho) * Xrt_Wr_Xr) + 1/sigma2_beta*I_p
  
    a <- y_1 - eta_f
    b <- (y_2[,2] - eta_r2) - rho * (y_2[,1]-eta_r1)
    mu_beta <- (1/sigma2 * (t(X_1)%*%W_1%*%a - 1/(1+rho) * t(X_2)%*%W_2%*%b))
    U_beta <- base::chol(Prec_beta)
    tmp_norm <- rnorm(nrow(Prec_beta))
    beta <- backsolve(U_beta, backsolve(U_beta, mu_beta, transpose=TRUE) + tmp_norm)
    
    ###sample eta_1
    ## need to group everyone who responded in week 1, so weeks_1 == 1 and weeks_2==2
    id1 <- which(weeks_1==1)
    id2 <- which(weeks_2==2)
    Psi_1f <- Psi_1[id1,]
    Psi_2r <- Psi_2[id2,]
    W_1f <- W_1[id1,id1]
    W_2r <- W_2[id2,id2]
    a <- y_1[id1] - X_1[id1,]%*%beta
    b <- -(y_2[id2,2] - X_2[id2,]%*%beta - Psi_2r %*% eta[,2] - rho*(y_2[id2,1] - X_2[id2,]%*%beta))
    Prec_eta1 <- 1/sigma^2 * (t(Psi_1f)%*%W_1f%*%Psi_1f + 
                       rho^2/(1-rho^2) * (t(Psi_2r)%*%W_2r %*% Psi_2r)) +  
                       (1/sigma2_eta_1 + phi^2/sigma2_eta) * I_D 
    mu_eta1 <- (1/sigma^2 * (t(Psi_1f) %*% W_1f %*% a + rho/(1-rho^2) * (t(Psi_2r)%*%W_2r%*% b)) +
                               phi/sigma2_eta*eta[,2])
    U_eta1 <- base::chol(Prec_eta1)
    tmp_norm <- rnorm(nrow(Prec_eta1))
    eta[,1] <- backsolve(U_eta1, backsolve(U_eta1, mu_eta1, transpose=TRUE) + tmp_norm)
    #sum to zero
    #eta[,1] <- eta[,1] - mean(eta[,1])
  
    ###sample all other etas
    for(t in 2:(n_weeks-1)){
      idf <- which(weeks_1==t) #indices for first-time at t
      idr <- which(weeks_2==t) #indices for followup at t
      idu <- which(weeks_2==(t+1)) #indices for followup at t+1

      Psi_tf <- Psi_1[idf,]
      Psi_tr <- Psi_2[idr,]
      Psi_tu <- Psi_2[idu,]

      W_tf <- W_1[idf,idf]
      W_tr <- W_2[idr,idr]
      W_tu <- W_2[idu,idu]

      a <- y_1[idf] - X_1[idf,]%*%beta
      b <- y_2[idr,2] - (1-rho) * X_2[idr,]%*%beta - rho*(y_2[idr,1] - Psi_tr%*%eta[,(t-1)])
      d <- Psi_tu%*%eta[,(t+1)] + (1-rho)*X_2[idu,]%*%beta - (y_2[idu,2] - rho*y_2[idu,1]) 
      
      Prec_eta_t <- (1+phi^2)/sigma2_eta * I_D + 
                         1/sigma^2 * (t(Psi_tf) %*% W_tf %*% Psi_tf +
                                      1/(1-rho^2)*(t(Psi_tr) %*% W_tr %*% Psi_tr +
                                                    rho^2 * (t(Psi_tu) %*% W_tu %*% Psi_tu )))
    
      mu_eta_t <-  (phi/sigma2_eta * (eta[,(t-1)] + eta[,(t+1)]) +
                                   1/sigma^2 * (t(Psi_tf) %*% W_tf %*% a +
                                                1/(1-rho^2) * (t(Psi_tr) %*% W_tr %*% b + rho * t(Psi_tu) %*% W_tu %*% d))
                                  )
    U_eta_t <- base::chol(Prec_eta_t)
    tmp_norm <- rnorm(nrow(Prec_eta_t))
    eta[,t] <- backsolve(U_eta_t, backsolve(U_eta_t, mu_eta_t, transpose=TRUE) + tmp_norm)
    }

    ###sample eta_T
    idf <- which(weeks_1==n_weeks)
    idr <- which(weeks_2==n_weeks)
    Psi_Tf <- Psi_1[idf,]
    Psi_Tr <- Psi_2[idr,]
    W_Tf <- W_1[idf,idf]
    W_Tr <- W_2[idr,idr]
    a <- y_1[idf] - X_1[idf,]%*%beta
    b <- y_2[idr,2] - X_2[idr,]%*%beta - rho*(y_2[idr,1] - X_2[idr,] %*% beta - Psi_Tr%*%eta[,(n_weeks-1)])
    
    Prec_eta_T <- 1/sigma2_eta * I_D + 1/sigma2 * (t(Psi_Tf) %*% W_Tf %*% Psi_Tf + 
                                                       1/(1-rho^2) * t(Psi_Tr) %*% W_Tr %*% Psi_Tr) 
    mu_eta_T <- (phi/sigma2_eta*eta[ ,(n_weeks-1)] + 1/sigma2 * (t(Psi_Tf) %*% W_Tf %*% a 
                                                                      + 1/(1-rho^2)*t(Psi_Tr) %*% W_Tr %*% b) )
  
    U_eta_T <- base::chol(Prec_eta_T)
    tmp_norm <- rnorm(nrow(Prec_eta_T))
    eta[,n_weeks] <- backsolve(U_eta_T, backsolve(U_eta_T, mu_eta_T, transpose=TRUE) + tmp_norm)
  
   #store this run
   if(iter > n_burn){#discard burn-in steps
     index <- iter-n_burn
     eta.post[index,,] <- as.matrix(eta)
     betas.post[index,] <- as.vector(beta)
     rho.post[index] <- rho
     phi.post[index] <- phi
     sigma2.post[index] <- sigma2
     sigma2_eta.post[index] <- sigma2_eta
     sigma2_eta_1.post[index] <- sigma2_eta_1
   }
  
  }
  posteriors <- list(eta=eta.post, betas=betas.post, rho=rho.post, phi=phi.post,
                     sigma2=sigma2.post, sigma2_eta=sigma2_eta.post,
                     sigma2_eta_1=sigma2_eta_1.post, rho.accept_rate=(rho.num_accept/n_keep))
  print(paste0("Rho accept rate = ", (rho.num_accept/n_keep)))
  return(posteriors)
}
