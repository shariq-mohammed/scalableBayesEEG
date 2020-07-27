#' Bayesian Estimation (Through Exact/Approximate Estimation) And Prediction
#'
#' Bayesian variable selection for the EEG data. The estimation type (exact or
#' approx) could needs to be specified. Also performs prediction for a test data
#' set. (Warning: type='exact' will require huge amount of time to execute even
#' for moderate values of \eqn{L} or \eqn{\tau}.)
#'
#' @param y binary response n values
#' @param X data - dimension \eqn{n}x\eqn{L}x\eqn{\tau}
#' @param c_init initial value for \eqn{c} - dimension \eqn{\tau} (defaults to 1)
#' @param beta_init initial values for \eqn{\beta} - dimension \eqn{L}x\eqn{\tau}
#'                  (defaults to 0)
#' @param zeta_init initial values for \eqn{\zeta} - dimension \eqn{L}x\eqn{\tau}
#'                  (defaults to randomly generated values from 1 or \eqn{v_0})
#' @param nusq_init initial values for \eqn{\nu^2} - dimension \eqn{L}x\eqn{\tau}
#'                  (defaults to 1)
#' @param l_init initial value for \eqn{\lambda} - dimension \eqn{L\tau}
#'               (defaults to 0)
#' @param prior.mu prior mean for \eqn{\lambda} - dimension \eqn{L\tau}
#' @param Spat.cov spatial covariance matrix  - dimension \eqn{L}x\eqn{L}
#' @param Temp.cov temporal covariance matrix  - dimension \eqn{\tau}x\eqn{\tau}
#' @param v0 hyperparameter \eqn{v_0}
#' @param a1 hyperparameter \eqn{a_1}
#' @param a2 hyperparameter \eqn{a_2}
#' @param q hyperparameter \eqn{q}
#' @param Nmcmc number of MCMC samples to generate
#' @param burnin number of samples to be dropped as burnin values
#' @param X.test Test data samples - dimension \eqn{m}x\eqn{L}x\eqn{\tau}
#' @param type Type of estimation (choose between 'exact' or 'approx')
#'
#' @return \itemize{
#' \item{\code{b.strucBayes}} {estimated coefficients \eqn{\beta} - dimension
#'                             \eqn{\tau}x\eqn{L}}
#' \item{\code{p.strucBayes.pred}} {estimated local prediction probabilities
#'                                  - dimension \eqn{m}x\eqn{\tau}}
#' \item{\code{y.strucBayes.pred}} {predicted response - length \eqn{m}}
#' \item{\code{run.time}} {total run time}
#' }
#'
#' @export
#' @importFrom stats rbinom kmeans median

spTempBayes = function(y, X, c_init = rep(1, tau), beta_init = NULL,
                       zeta_init = NULL, nusq_init = NULL, l_init = NULL,
                       prior.mu, Spat.cov, Temp.cov, v0, a1, a2, q,
                       Nmcmc, burnin, X.test = NULL, type){
  n = dim(X)[1]
  p = dim(X)[2]
  tau = dim(X)[3]
  m = dim(X.test)[1]

  # scale data at subject level by its Frobenius norm
  X_int = sapply(1:n, function(i) sqrt(sum(X[i,,]^2)))
  X_sc = array(NA, dim = dim(X))
  for(i in 1:n) X_sc[i,,] = X[i,,]/X_int[i]

  if(!is.null(m)){
    X.test_int = sapply(1:m, function(i) sqrt(sum(X.test[i,,]^2)))
    X.test_sc = array(NA, dim = dim(X.test))
    for(i in 1:m) X.test_sc[i,,] = X.test[i,,]/X.test_int[i]
  }

  if(is.null(beta_init)) beta_init = matrix(0, nrow = p, ncol = tau)
  if(is.null(zeta_init)){
    zeta_init = matrix(rbinom(p*tau,size = 1, prob = 0.5), nrow = p, ncol = tau)
    zeta_init[zeta_init==0] = v0
  }
  if(is.null(nusq_init)) nusq_init = matrix(1, nrow = p, ncol = tau)
  if(is.null(l_init)) l_init = rep(0, p*tau)

  ind = seq(burnin+1, Nmcmc, by=10)

  if(type == 'exact'){
    stm = proc.time()
    sim = spTempExact(y, X_sc, c_init, beta_init, zeta_init, nusq_init,
                      l_init, prior.mu, Spat.cov, Temp.cov, v0, a1, a2,
                      q, Nmcmc, ind)
    run.time = proc.time()-stm
  }else if(type == 'approx'){
    stm = proc.time()
    sim = spTempApprox(y, X_sc, c_init, beta_init, zeta_init, nusq_init,
                       l_init, prior.mu, Spat.cov, Temp.cov, v0, a1, a2,
                       q, Nmcmc, ind)
    run.time = proc.time()-stm
  } else{
    stm = proc.time()
    print('Specify the type of estimation procedure: approx or exact!')
    run.time = proc.time()-stm
  }

  loc.probs = matrix(nrow = tau, ncol = p)
  for(t in 1:tau) loc.probs[t,] = colMeans(sim$zeta[,,t]==1)
  sort.probs = apply(loc.probs, 2, sort)

  sort.clus = kmeans(t(sort.probs), centers = 2)
  ind1 = which(sort.clus$cluster==1)
  ind2 = which(sort.clus$cluster==2)
  clus1.mean = mean(loc.probs[,ind1])
  clus2.mean = mean(loc.probs[,ind2])

  b.strucBayes = matrix(0, nrow = tau, ncol = p)
  if(clus1.mean==clus2.mean) b.strucBayes = apply(sim$b, c(3,2), mean)
  if(clus1.mean>clus2.mean) b.strucBayes[,ind1] = apply(sim$b,c(3,2),mean)[,ind1]
  if(clus2.mean>clus1.mean) b.strucBayes[,ind2] = apply(sim$b,c(3,2),mean)[,ind2]

  c.est = apply(sim$c, 2, median)

  p.strucBayes.pred = NULL
  y.strucBayes.pred = NULL

  if(!is.null(X.test)){
    p.strucBayes.pred = pnorm(sapply(1:tau,
                                     function(t){
                                       crossprod(t(X.test_sc[,,t]),
                                                 b.strucBayes[t,])/c.est[t]
                                     }))
    y.matrix.pred = matrix(nrow = m, ncol = tau)
    for(i in 1:m) y.matrix.pred[i,] = as.numeric(p.strucBayes.pred[i,]>0.5)

    w = matrix(nrow = m, ncol = tau)
    for(i in 1:m){
      for(t in 1:tau){
        if(p.strucBayes.pred[i,t]>0.5) w[i,t] = p.strucBayes.pred[i,t]-0.5
        if(p.strucBayes.pred[i,t]<0.5) w[i,t] = 0.5-p.strucBayes.pred[i,t]
      }
      w[i,] = (w[i,]^2)/sum(w[i,]^2)
    }

    y.w = sapply(1:m, function(i) sum(w[i,]*y.matrix.pred[i,]))
    y.strucBayes.pred = as.numeric(y.w>0.5)
  }

  list(b.strucBayes = b.strucBayes,
       p.strucBayes.pred = p.strucBayes.pred,
       y.strucBayes.pred = y.strucBayes.pred,
       run.time = run.time)
}
