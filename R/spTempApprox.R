#' MCMC Algorithm For Estimation Of Local Model Through Scalable (Approximate)
#' Estimation
#'
#' Gibbs sampling, and Metropolis-Hastings within Gibbs sampling for a local
#' model through scalable (approximate) estimation.
#'
#' @param y binary response n values
#' @param X data - dimension \eqn{n}x\eqn{L}x\eqn{\tau}
#' @param c_init initial value for \eqn{c} - dimension \eqn{\tau}
#' @param beta_init initial values for \eqn{\beta} - dimension \eqn{L}x\eqn{\tau}
#' @param zeta_init initial values for \eqn{\zeta} - dimension \eqn{L}x\eqn{\tau}
#' @param nusq_init initial values for \eqn{\nu^2} - dimension \eqn{L}x\eqn{\tau}
#' @param l_init initial value for \eqn{\lambda} - dimension \eqn{L\tau}
#' @param prior.mu prior mean for \eqn{\lambda} - dimension \eqn{L\tau}
#' @param Spat.cov spatial covariance matrix  - dimension \eqn{L}x\eqn{L}
#' @param Temp.cov temporal covariance matrix  - dimension \eqn{\tau}x\eqn{\tau}
#' @param v0 hyperparameter \eqn{v_0}
#' @param a1 hyperparameter \eqn{a_1}
#' @param a2 hyperparameter \eqn{a_2}
#' @param q hyperparameter \eqn{q}
#' @param Nmcmc number of MCMC samples to generate
#' @param ind indices of MCMC samples to use
#'
#' @return \itemize{
#' \item{\code{b}} {MCMC samples for \eqn{\beta} - dimension (size of
#'                  \code{ind})x\eqn{L}x\eqn{\tau}}
#' \item{\code{zeta}} {MCMC samples for \eqn{\zeta} - dimension (size of
#'                  \code{ind})x\eqn{p}x\eqn{\tau}}
#' \item{\code{l}} {MCMC samples for \eqn{\lambda} - dimension (size of
#'                  \code{ind})x\eqn{L\tau}}
#' \item{\code{c}} {MCMC samples for \eqn{c} - dimension (size of
#'                  \code{ind})x\eqn{L\tau}}
#' \item{\code{t0}} {Samples for \eqn{t0} - dimension (size of
#'                  \code{ind})}
#' }
#'
#' @export
#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rbinom rgamma

########################################################
### MCMC algorithm for the estimation of local model ###
########################################################

spTempApprox = function(y, X, c_init, beta_init, zeta_init, nusq_init,
                        l_init, prior.mu, Spat.cov, Temp.cov,
                        v0, a1, a2, q, Nmcmc, ind ){
  n = dim(X)[1]
  p = dim(X)[2]
  tau = dim(X)[3]

  c = c_init
  b = beta_init
  zeta = zeta_init
  nusq = nusq_init
  l = l_init

  # Storing all the required parameters from the MCMC samples
  c_store = matrix(nrow = Nmcmc, ncol = tau)
  zeta_store = array(dim = c(Nmcmc, p, tau))
  b_store = array(dim = c(Nmcmc, p, tau))
  l_store = matrix(nrow = Nmcmc, ncol = p*tau)
  t0_store = numeric(length = Nmcmc)

  Spat.inv = chol2inv(chol(Spat.cov))

  # posterior mean of latent parameter z
  zMean = sapply(1:tau, function(t) crossprod(t(X[,,t]), b[,t]))

  # start Gibbs sampling
  for(i in 1:Nmcmc){
    if(i %% 100 == 0) print(paste("Running iteration # = ", i, sep=""))
    z = sapply(1:tau,
               function(t){
                 (y*rtruncnorm(1, 0, Inf, zMean[,t], sd = sqrt(c[t]))+
                    (1-y)*rtruncnorm(1, -Inf, 0, zMean[,t], sd = sqrt(c[t])))
               })
    # update beta
    b = sapply(1:tau,
               function(t){
                 symX = crossprod(X[,,t], X[,,t]/c[t])
                 gm = nusq[,t]*zeta[,t]
                 GMinv = diag(1/gm)
                 bVar = chol2inv(chol(symX+GMinv))
                 temp1 = crossprod(X[,,t], z[,t]/c[t])
                 bMean = crossprod(t(bVar), temp1)

                 rmvnorm(1, mean = bMean, sigma = bVar)
               })
    b_store[i,,] = b

    # update zeta
    zeta = sapply(1:tau,
                  function(t){
                    tempExp = exp(-(b[,t]^2)/(2*nusq[,t]))
                    pr1 = (1-pnorm(l[((t-1)*p)+(1:p)]))*(tempExp^(1/v0))/sqrt(v0)
                    pr2 = pnorm(l[((t-1)*p)+(1:p)])*tempExp

                    temp.probs = pr2/(pr1+pr2)
                    temp.probs[is.na(temp.probs)] = 0.5
                    zet = rbinom(p, 1, prob = temp.probs)
                    zet[zet==0] = v0
                    zet
                  })
    zeta_store[i,,] = zeta

    t0 = sample.int(tau, size = 1)
    t0_store[i] = t0
    # update lambda using M-H within Gibbs sampling
    for(k in 1:p){
      temp.ind = p*(t0-1)
      cond.var = Temp.cov[t0,t0]/Spat.inv[k,k]
      cond.mu = prior.mu[temp.ind+k]+crossprod(-Spat.inv[k,-k]/Spat.inv[k,k],
                                               (l[(1:p)+temp.ind]-
                                                  prior.mu[(1:p)+temp.ind])[-k])

      if(zeta[k,t0]==1) l[temp.ind+k] = metropolis(l[temp.ind+k], cond.mu,
                                                   sqrt(cond.var), 1)
      if(zeta[k,t0]==v0) l[temp.ind+k] = metropolis(l[temp.ind+k], cond.mu,
                                                    sqrt(cond.var), -1)
    }
    l.ind = (1:(p*tau))[-((p*(t0-1))+(1:p))]
    l.temp = sapply(l.ind,
                    function(lt){
                      t = ceiling(lt/p)
                      k = lt - (p*(t-1))
                      temp.ind = p*(t0-1)
                      cond.mu = prior.mu[lt]+((Temp.cov[t,t0]/Temp.cov[t0,t0])*
                                                (l[k+temp.ind]-
                                                   prior.mu[k+temp.ind]))
                      cond.var = (Temp.cov[t,t]-(Temp.cov[t,t0]*Temp.cov[t0,t]/
                                                   Temp.cov[t0,t0]))*
                        Spat.cov[k,k]
                      if(zeta[k,t]==1) return(metropolis(l[lt], cond.mu,
                                                         sqrt(cond.var), 1))
                      if(zeta[k,t]==v0) return(metropolis(l[lt], cond.mu,
                                                          sqrt(cond.var), -1))
                    })
    l[l.ind] = l.temp
    l_store[i,] = l

    # update nu^2
    nusq = sapply(1:tau,
                  function(t) 1/rgamma(p, shape = a1+0.5,
                                       rate = a2+((b[,t]^2)/(2*zeta[,t]))))

    # update c
    zMean = sapply(1:tau, function(t) crossprod(t(X[,,t]), b[,t]))

    c = sapply(1:tau, function(t) 1/rgamma(1, shape = (n+q)/2,
                                           rate = (q+sum((z[,t]-zMean[,t])^2))/2))
    c_store[i,] = c
  }

  # return samples of beta, zeta, lambda and c as they will be needed for ...
  # ... variable selection and prediction
  list(b = b_store[ind,,],
       zeta = zeta_store[ind,,],
       l = l_store[ind,],
       c = c_store[ind,],
       t0 = t0_store[ind])
}
