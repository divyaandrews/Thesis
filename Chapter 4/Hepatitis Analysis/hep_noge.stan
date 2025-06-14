functions {
   // Custom function to calculate the indicator function indf
    int indf(int x) {
        return (x == 0) ? 1 : 0;
    }

  real novel_geometric_lpmf(int x, real phi, real theta) {
    // Generalized Poisson log-probability mass function
    return indf(x) * log(phi) +
                   (1 - indf(x)) * (log(1 - phi) + (x - 1) * log(1 - theta) + log(theta));
  }
}

data {
  int<lower=0> N; // Number of time periods
  int Y[N];   // The time series data
}

parameters {
  real<lower=0> alpha0;          // alpha0
  real beta1_unconstd;           // Unconstrained beta1
  real alpha1_unconstd;          // Unconstrained alpha1
  real phi_unconstd;            // Unconstrained phi
  //real<lower=0> n0;             // n0
  real<lower=0> m1;              // Initial value for cmean
}

transformed parameters {
  real alpha1 = inv_logit(alpha1_unconstd);  // Transformed alpha1
  real beta1 = 0.62 * inv_logit(beta1_unconstd); // Transformed beta1
  real phi = inv_logit(phi_unconstd);  // Transformed phi
  real cmean[N];
    // Initialize the conditional mean
  cmean[1] = m1; // Initial value for cmean

  for (n in 2:N) {
    cmean[n] = alpha0 + alpha1 * Y[n-1] + beta1 * cmean[n-1];
  }
}

model {
  // Priors for alpha0
  alpha0 ~ lognormal(0.4131, 0.3);
  // Priors for theta
  phi_unconstd ~ normal(0.9639, 50);
  
  // Bivariate normal prior on transformed unconstrained parameters
   vector[2] mu = [ -0.487, -0.9048]';           // Mean vector
  matrix[2, 2] Sigma = [[100, 50],  // Covariance matrix
                        [50, 100]];
  
  // Apply bivariate normal prior
  [alpha1_unconstd, beta1_unconstd] ~ multi_normal(mu, Sigma);
  // Prior for m1
  m1 ~ lognormal(2.385, 0.3);
  
  // Likelihood
  for (n in 2:N) {
    Y[n] ~ novel_geometric(phi, (1-phi)/cmean[n]);  // Likelihood for Generalized Poisson
  }
}

generated quantities {
  real log_lik; 
  real cmean_last;// Sum of log-likelihood
  log_lik = 0;
  for (n in 1:N) {
    log_lik += novel_geometric_lpmf(Y[n] | phi, (1-phi)/cmean[n]);  // Summing Poisson log-likelihoods
  }
  cmean_last = cmean[N];
}
