###############################################################################
# log prob of full conditionals of parameters
###############################################################################

log_p_M_t = function(M_t, rho, y_t, x_c_t){
  # M_t: proposed value for M_t (real)
  # rho: current default correlation (real in (0,1))
  # y_t: y for time t (vector length N)
  # x_c_t: critical values for time t (vector length N)
  # output: log prob full conditional for M_t

  phi_e_star = pnorm((x_c_t-M_t*sqrt(rho))/(sqrt(1-rho)))

  lp=sum(ifelse(y_t == 1L, log(phi_e_star), log(1-phi_e_star)))
  lp=lp-0.5*M_t*M_t
  return(lp)
}

log_p_rho = function(rho, M, y, x_c){
  # rho: current default correlation (real in (0,1))
  # M: common factors (length Tau)
  # y: observations (matrix size N\times Tau ints in {0,1})
  # x_c: critical values (matrix size N\times Tau)
  # output: log prob full conditional for rho

  phi_e_star = sweep(x=x_c,MARGIN=2L,STATS=M*sqrt(rho))
  phi_e_star = pnorm(phi_e_star/sqrt(1-rho))

  lp=sum(ifelse(y == 1L, log(phi_e_star), log(1-phi_e_star)))
  return(lp)
}

###############################################################################
# MH within Gibbs updates
###############################################################################

# function for getting shape1 for dbeta based on shape2 and mode
get_a = function(m, b) {(m*(b-2)+1)/(1-m)}

update_rho = function(rho, M, y, x_c, b){
  # rho: current default correlation (real in (0,1))
  # M: common factors (length Tau)
  # y: observations (matrix size N\times Tau ints in {0,1})
  # x_c: critical values (matrix size N\times Tau)
  # b: parameter shape2 for the proposal (higher -> more concentrated on rho)
  # output: new rho

  # proposal distribution sets mode at current value
  r = rbeta(n=1L, shape1 = get_a(rho,b), shape2 = b) # proposal
  lp = log_p_rho(rho=r, M=M, y=y, x_c=x_c) # full conditional at proposal
  lq = dbeta(r, shape1 = get_a(rho,b), shape2 = b, log = TRUE) # density at proposal
  lp_stay = log_p_rho(rho=rho, M=M, y=y, x_c=x_c) # full conditional at stay
  lq_stay = dbeta(rho, shape1=get_a(r,b), shape2=b, log=TRUE) # density of reverse proposal
  l_a_ratio = min(0, (lp-lq)-(lp_stay-lq_stay))
  p_accept = exp(l_a_ratio)
  accept = (runif(1L) <= p_accept)
  accept = ifelse(is.na(accept), FALSE, accept) # if p_accept is NaN, accept is NA

  if (accept) {
    r
  } else{
    rho
  }

}

update_M_t = function(M_t, rho, y_t, x_c_t){
  # M_t: current value for M_t (real)
  # rho: current default correlation (real in (0,1))
  # y_t: y for time t (vector length N)
  # x_c_t: critical values for time t (vector length N)
  # output: new M_t

  m = rnorm(n = 1L, mean = M_t) # proposal
  lp = log_p_M_t(M_t=m, rho=rho, y_t=y_t, x_c_t=x_c_t)
  lp_stay = log_p_M_t(M_t=M_t, rho=rho, y_t=y_t, x_c_t=x_c_t)
  p_accept = min(1, exp(lp-lp_stay)) # proposal kernel is reversible
  accept = (runif(1L) <= p_accept)
  accept = ifelse(is.na(accept), FALSE, accept) # if p_accept is NaN, accept is NA

  if(accept){
    m
  }else{
    M_t
  }

}

###############################################################################
# iterations
###############################################################################

# initialize chain
initialize_chain = function(Tau){
  # Tau: number of timesteps
  # output: vector of initial points

  # initalize latent variables
  values = numeric(Tau+1L)
  values[1L:Tau] = rnorm(Tau)
  values[Tau+1L] = rbeta(n=1L, shape1 = get_a(0.5,10), shape2 = 10)

  return(values)
}

#' Bayesian Inference of Default Correlations Using Probabilities of Default
#'
#' This functions runs an MCMC chain to perform Bayesian inference on the default
#' correlation parameter, by using probabilites of default estimated for each
#' client. See Details for more information.
#'
#' The MCMC chain is a Metropolis-Hastings-withing-Gibbs (MHwG) sampler. The MH
#' updates are needed because the full conditionals of the parameters are
#' not analytically tractable. The prior on rho is uniform in (0,1), and for the
#' latent factors are standard normals.
#'
#' The approach used here requires knowing the actual outcomes of the defaults
#' for a fixed number of clients (N) over a time period (tau), which are encoded
#' in the matrix \code{y}. It also requires an estimate of the probability of
#' default (PD) for each entry in \code{y} in a matrix we call \code{p_d}.
#' The inferences are better when the estimates of the PDs are better.
#'
#' It is crucial to understand that client n in time t1 does not need to
#' correspond to the same person in the n-th row at time t2. In other words,
#' we assume exchangeability within the rows of y. The user only needs to
#' make sure that every entry in \code{p_d} is an estimate for the corresponding
#' entry in \code{y}.
#'
#' @param S number of iterations for the chain.
#' @param y integer matrix of size \code{N x tau} containing a 1 in (n,t)
#' if the n-th client at time t defaulted, and 0 otherwise.
#' @param p_d matrix with probability of default of the n-th client at
#' time t.
#' @param b shape2 parameter of the beta proposal for rho. Default value is 50
#' tends to work well.
#' @param verbose integer. Function prints every \code{verbose} iterations. The
#' default value of 0L prints no output.
#' @param init an optional initialization for the chain. If \code{NULL}
#' (default), the M_t's are initialized from the prior and rho from a beta close
#' to 0.5.
#'
#' @return a matrix of size \code{S x (tau+1)}, where the first tau columns
#' correspond to the latent factors, and the last one contains samples for rho.
#'
#' @export
bidc_pd = function(S = 2000L, y, p_d, b = 50, verbose = 0L, init = NULL){

  Tau = ncol(y) # number of timesteps
  stopifnot(0 <= min(p_d) && max(p_d) <= 1) # sanity check
  x_c = qnorm(p_d) # inverse probit transform
  # y = y + 1L # use internal representation of (1,2)
  stopifnot(Tau == ncol(x_c))

  # define matrix of values of chain
  chain = matrix(NA_real_, nrow = S, ncol = Tau+1L)

  # initialize
  if(is.null(init)){
    chain[1L,] = initialize_chain(Tau=Tau)
  } else{
    chain[1L,] = init
  }

  # loop
  for(s in 2L:S){
    # update M
    for(tau in 1L:Tau){
      chain[s, tau] = update_M_t(
        M_t   = chain[s - 1L, tau],
        rho   = chain[s - 1L, Tau + 1L],
        y_t   = y[, tau],
        x_c_t = x_c[, tau]
      )
    }

    # update rho
    chain[s, Tau + 1L] = update_rho(
      rho    = chain[s - 1L, Tau + 1L],
      M      = chain[s, 1L:Tau],
      y      = y,
      x_c    = x_c,
      b      = b
    )

    # print status
    if(verbose>0L && (s %% verbose == 0L)){
      cat(sprintf("Iteration %d: rho=%.2f, M_1=%.2f, M_%d=%.2f\n",
                  s, chain[s, Tau + 1L], chain[s, 1L], Tau, chain[s, Tau]))
    }

  }

  return(chain)

}
