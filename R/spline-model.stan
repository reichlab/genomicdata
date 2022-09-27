//
// Stan program for fitting splines to compositional count data of number of 
// samples per clade over time.
// 
//

data {
  int<lower=1> T; // # of observations
  int<lower=1> K; // # of clades
  
  array[T,K] int Y; // counts of clade k at time t
  array[T] int N;   // counts of total samples at each time t
  
  int<lower=1> n_basis; // number of basis functions
  matrix[T, n_basis] B; // design matrix with basis functions **will increase T when forecasting**
  
  real<lower=0> sigma_beta; // prior sd for beta matrix params
  real<lower=0> sigma_beta_rw; // prior scale for beta cauchy RW
  real<lower=0> A_eps; // prior parameter for clade-specific sd for noise around phi, assuming Gamma(2, 1/A), based on Chung (2013)
}

transformed data {
  vector[K-1] ones;
  for (i in 1:(K-1)){
    ones[i] = 1;
  }
}

parameters {
  matrix[n_basis,K-1] beta_raw;
  real<lower=0> sigma_eps;
  matrix[T,K] eps;
}

transformed parameters{
  matrix[T,K] phi; 
  array[T] simplex[K] pi;
  matrix[n_basis,K] beta = append_col(beta_raw, -beta_raw*ones);

  // compute spline values, for each k (clade)
  for (k in 1:K){
    phi[,k] = B * beta[,k] + eps[,k];
  }
  
  // compute probabilities, for each time
  for (t in 1:T){
    pi[t] = softmax(to_vector(phi[t,]));
    //eps[t,] = to_row_vector(pi[t] - to_vector(Y[t,])/N[t]);
  }
  
  // for generating posterior phi
  // matrix[T,K] phi_mean; 
  // matrix[T,K] phi_lb; 
  // matrix[T,K] phi_ub; 
  // phi_mean = phi - eps;
  // for (k in 1:K){
  //   phi_lb[,k] = phi_mean[,k] - 1.28*sigma_eps[k];
  //   phi_ub[,k] = phi_mean[,k] + 1.28*sigma_eps[k];
  // }
}

model {
  // RW-like priors
  to_vector(beta_raw[1,1:(K-1)]) ~ normal(0, sigma_beta); // start at zero
  for (j in 2:n_basis){
    for (k in 2:(K-1)){
      beta[j,k] ~ cauchy(beta[j-1,k], sigma_beta_rw); // RW-like structure after 
    }
  }
  
  // distribution of random noise
  sigma_eps ~ gamma(2, 1/A_eps);
  to_vector(eps) ~ normal(0, sigma_eps);
  
  // define probabilities and define likelihood, for each time
  for (t in 1:T){
    Y[t,] ~ multinomial(pi[t]);
  }
}

generated quantities {
  // for simulating new observations
  // array[T,K] int y_sim;
  // for (t in 1:T){
  //   y_sim[t,] = multinomial_rng(pi[t], N[t]);
  // }
}

