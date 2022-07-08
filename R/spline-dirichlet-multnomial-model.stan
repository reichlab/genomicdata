//
// Stan program for fitting splines to compositional count data of number of 
// samples per clade over time.
// 
//

data {
  int<lower=1> T; // # of observations
  int<lower=1> K; // # of clades
  
  array[T,K] int Y; // counts of clade k at time t
  
  int<lower=1> n_basis; // number of basis functions
  matrix[T, n_basis] B; // design matrix with basis functions **will increase T when forecasting**
  
  real<lower=0> sigma_beta; // prior sd for beta matrix params
  real<lower=0> sigma_beta_rw; // prior scale for beta cauchy RW
}

transformed data {
  // vector used to define beta_K as the sum of all other betas
  vector[K-1] ones;
  for (i in 1:(K-1)){
    ones[i] = 1;
  }

  matrix<lower=0,upper=1>[T,K] F; // frequencies of clade k at time t
  for (t in 1:T) {
    F[t,] = to_row_vector(Y[t,]) / sum(Y[t,]);
  }
}

parameters {
  matrix[n_basis,K-1] beta_raw;
}

transformed parameters{
  matrix[T,K] phi; 
  
  // defines beta as beta_raw plus a column that 
  matrix[n_basis,K] beta = append_col(beta_raw, -beta_raw*ones);

  // compute spline values, for each k (clade)
  for (k in 1:K){
    phi[,k] = exp( B * beta[,k] );
  }
}

model {
  // RW-like priors
  // to_vector(beta_raw[1,1:(K-1)]) ~ normal(0, sigma_beta); // start at zero
  // for (j in 2:n_basis){
  //   for (k in 2:(K-1)){
  //     beta[j,k] ~ cauchy(beta[j-1,k], sigma_beta_rw); // RW-like structure after 
  //   }
  // }
  to_vector(beta_raw[,1:(K-1)]) ~ normal(0, sigma_beta);
  
  // define probabilities and define likelihood, for each time
  for (t in 1:T){
    F[t,] ~ dirichlet(phi[t,]);
  }
}

generated quantities {
  // for simulating new observations
  matrix<lower=0,upper=1>[T,K] F_sim;
  for (t in 1:T){
    F_sim[t,] = to_row_vector(dirichlet_rng(to_vector(phi[t,])));
  }
}

