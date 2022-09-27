// Stan program for simple dirichlet "model 4" - with pis for each obs

data {
  int<lower=1> T; // # of observations
  int<lower=1> K; // # of clades
  
  array[T,K] int Y; // frequencies of clade k at time t

  int<lower=1> n_basis; // number of basis functions
  matrix[T, n_basis] B; // design matrix with basis functions

  real<lower=0> mu_beta;    // prior mean for beta
  real<lower=0> sigma_beta; // prior sd for beta
}

transformed data {
  // vector used to define beta_K as the sum of all other betas
  vector[K-1] ones;
  for (i in 1:(K-1)){
    ones[i] = 1;
  }
}

parameters {
  matrix<lower=0>[n_basis,K] beta; // note truncation >0
  array[T] simplex[K] pi;
}

transformed parameters{
  // compute spline values: can't be negative due to basis>0 and beta>0 above
  matrix<lower=0>[T,K] phi = B * beta ;
}

model {
  // priors
  to_vector(beta) ~ cauchy(mu_beta, sigma_beta);

  // define likelihood, for each time
  for (t in 1:T){
    pi[t] ~ dirichlet(phi[t,]);
    Y[t,] ~ multinomial(pi[t]);
  }
}

generated quantities {
  // for simulating new observations
  matrix<lower=0,upper=1>[T,K] y_sim;
  for (t in 1:T){
    y_sim[t,] = to_row_vector(dirichlet_rng(to_vector(pi[t,])));
  }
}

