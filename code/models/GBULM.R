library(mvtnorm)

model1 <- function(y, X, Psi, wgt, iter=500, burn=200){
  n <- length(y)
  p <- ncol(X)
  r <- ncol(Psi)
  WGT <- Diagonal(n, wgt)
  XWX <- t(cbind(X, Psi))%*%WGT%*%cbind(X,Psi)
  XWY <- t(cbind(X, Psi))%*%WGT%*%y
  sig2<- 1
  sig2RE <- 1
  sig2Out <- rep(NA, iter)
  betaOut <- matrix(NA, nrow=iter, ncol=p+r)
  #pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    if(i %% 100 == 0){print(paste("Gibbs indep. iter", i))}
    A <- XWX/sig2 + bdiag(Diagonal(p, 1/1000), Diagonal(r, 1/sig2RE))
    Ainv <- solve(A)
    b <- XWY/sig2
    beta <- betaOut[i,] <- as.vector(rmvnorm(1,Ainv%*%b,as.matrix(Ainv)))

    sig2 <- sig2Out[i] <- 1/rgamma(1,
                                   .1 + n/2,
                                   .1 + 0.5*crossprod(wgt*(y-cbind(X, Psi)%*%beta), (y-cbind(X, Psi)%*%beta)))
    sig2RE <- 1/rgamma(1,
                       .1 + r/2,
                       .1 + 0.5*sum(beta[-c(1:p)]^2))
    #setTxtProgressBar(pb, i)
  }
  return(list(Beta=betaOut[-c(1:burn),], sig2=sig2Out[-c(1:burn)]))
}

fit_independent_gibbs <- function(){
    ests <- c()
    for(t in 1:n_weeks){
      HPS_week <- HPS_sample_long %>% filter(WEEK==t)
      X_week <- model.matrix(~COVAR.1 + COVAR.2 + I(COVAR.2^2), data=HPS_week)
      N_week <- nrow(HPS_week)
      scale_weights <- HPS_week$PWEIGHT * N_week/sum(HPS_week$PWEIGHT)

      mod <- model1(y=HPS_week$RESPONSE, X=X_week, 
                    Psi=model.matrix(~AREA-1, data=HPS_week), 
                    wgt=scale_weights, iter=iter, burn=burn)
      popWeek <- HPS_pop_long %>% filter(WEEK==t)
      popX <- model.matrix(~COVAR.1 +COVAR.2+I(COVAR.2^2), data=popWeek)
      popPsi <- model.matrix(~AREA-1, data=popWeek)
      preds <- cbind(popX, popPsi) %*% t(mod$Beta)


      sigs <- matrix(rep(sqrt(mod$sig2), nrow(popX)), nrow=nrow(popX), byrow=T)

      preds <- rnorm(n=length(preds),
                     mean=c(preds),
                     sd=c(sigs))
      preds <- matrix(preds,
                      nrow=nrow(popX))

      agg <- data.frame(AREA=popWeek$AREA, preds) %>%
        group_by(AREA) %>%
        summarize_all(mean)
      ests <- c(ests, rowMeans(agg[,-1]))
      rm(popX, popPsi, popWeek)
    }
}
