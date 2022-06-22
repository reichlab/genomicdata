## simulate clade data from Stan spline model

library(tidyverse)
library(cmdstanr)
library(posterior)
library(splines)

mdl <- cmdstan_model(stan_file = "R/spline-model.stan")

########
## simulate fake data
T <- 100
K <- 3
t <- 1:T
n_basis <- round(T/10)
B <- bs(t, df = n_basis, intercept=TRUE)
N <- rep(200, times=T)  

# simulate betas
sigma_beta <- 1
sigma_beta_rw <- 0.2
sim_betas <- function(n_basis, K, s_beta, s_beta_rw) {
  beta <- matrix(NA, nrow=n_basis, ncol=K)
  beta[1,] <- rnorm(K, 0, s_beta)
  for(j in 2:n_basis) {
    beta[j,] <- beta[j-1, ] + rcauchy(K, 0, s_beta_rw)
  }
  return(beta)
}
nsim <- 20
set.seed(731)
betas <- replicate(nsim, 
                       sim_betas(n_basis=n_basis, K=K, s_beta=sigma_beta, s_beta_rw=sigma_beta_rw), 
                       simplify=FALSE)

# simulate meaningless Ys
weights <- B %*% beta_list[[1]]
Y <- matrix(NA, nrow=T, ncol=K)
for(t in 1:T)
  Y[t,] <- t(rmultinom(1, size=N, prob=exp(weights[t,])/sum(exp(weights[t,])))) #matrix(10, nrow=T, ncol=K)

sim_out <- mdl$sample(data = list(T = T, ## num observations
                                K = K,   ## num clades
                                # Y should be real but fake data if we're simulating
                                # new Y
                                Y = Y,    ## obs
                                N = N,    ## total samples per time
                                n_basis = n_basis, ## num basis funcitons
                                B = B,        ## design matrix of basis functions 
                                sigma_beta = 5, ## prior sd for beta
                                sigma_beta_rw = 5), ## prior sd for beta RW 
                      init = list(
                        list(beta = betas[[1]]),list(beta = betas[[2]]),
                        list(beta = betas[[3]]),list(beta = betas[[4]]),
                        list(beta = betas[[5]]),list(beta = betas[[6]]),
                        list(beta = betas[[7]]),list(beta = betas[[8]]),
                        list(beta = betas[[9]]),list(beta = betas[[10]]),
                        list(beta = betas[[11]]),list(beta = betas[[12]]),
                        list(beta = betas[[13]]),list(beta = betas[[14]]),
                        list(beta = betas[[15]]),list(beta = betas[[16]]),
                        list(beta = betas[[17]]),list(beta = betas[[18]]),
                        list(beta = betas[[19]]),list(beta = betas[[20]])
                      ),
                      fixed_param = TRUE,
                      chains=20, 
                      parallel_chains=10)

## get and plot simulated ys
draws_y <- as_draws_df(sim_out$draws("y_sim"))

ysim <- pivot_longer(draws_y, cols = starts_with("y_sim")) %>%
  rename(chain=.chain) %>%
  group_by(name, chain) %>%
  summarize(med = median(value),
            q10 = quantile(value, probs=0.1),
            q90 = quantile(value, probs=0.9)) %>%
  tidyr::extract(name, c("week", "clade"), "([[:alnum:]]+),([[:alnum:]]+)", remove=FALSE, convert=TRUE)
  
ggplot(ysim, aes(x=week, y=med)) +
  geom_line(aes(color=factor(clade))) +
  facet_wrap(.~chain) +
  ylab("clade counts")


## run HMC for one beta and compare to truth
idx <- 13
Y <- ysim %>% 
  filter(chain==idx) %>%
  arrange(week) %>%
  pivot_wider(id_cols = week, names_from = "clade", values_from = "med", names_prefix = "clade") %>%
  as.matrix() %>%
  .[,2:(K+1)]
  
mdl_fit <- mdl$sample(data = list(T = T, ## num observations
                                  K = K,   ## num clades
                                  Y = Y,    ## obs
                                  N = N,    ## total samples per time
                                  n_basis = n_basis, ## num basis funcitons
                                  B = B,        ## design matrix of basis functions 
                                  sigma_beta = 5, ## prior sd for beta
                                  sigma_beta_rw = 2), ## prior sd for beta RW 
                      chains=4, 
                      parallel_chains=4)

## get and plot simulated pi's
draws_pi <- as_draws_df(mdl_fit$draws("pi"))

pi_est <- pivot_longer(draws_pi, cols = starts_with("pi")) %>%
  rename(chain=.chain) %>%
  group_by(name) %>%
  summarize(med = median(value),
            q10 = quantile(value, probs=0.1),
            q90 = quantile(value, probs=0.9)) %>%
  tidyr::extract(name, c("week", "clade"), "([[:alnum:]]+),([[:alnum:]]+)", remove=FALSE, convert=TRUE)

true_pis <- B %*% betas[[idx]] %>%
  as_tibble() %>%
  mutate(pi_1 = exp(V1) / ( exp(V1)+exp(V2)+exp(V3) ),
         pi_2 = exp(V2) / ( exp(V1)+exp(V2)+exp(V3) ),
         pi_3 = exp(V3) / ( exp(V1)+exp(V2)+exp(V3) )) %>%
  select(starts_with("pi")) %>%
  mutate(week = 1:n()) %>%
  pivot_longer(cols = starts_with("pi"), names_to = "clade", values_to = "pi")


ggplot(pi_est, aes(x=week)) +
  geom_line(aes(y=med, color=factor(clade))) + 
  geom_line(data=true_pis, aes(y=pi, group=clade)) +
  geom_ribbon(aes(fill=factor(clade), ymin=q10, ymax=q90), alpha=.3)+
  ylim(c(0,1)) +
  ylab("clade probabilities")


########
## estimate model using real data
clades <- c("20A", "20C", 
            "20G", "20I (Alpha, V1)", 
            "21I (Delta)", "21J (Delta)", 
            "21K (Omicron)")
raw_d <- read_csv("data/genomic_data_for_modeling.csv") %>%
  filter(Nextstrain_clade %in% clades, 
         epidate <= as.Date("2021-06-01")) 

d <- raw_d %>%
  pivot_wider(id_cols = epiweek_year, 
              names_from = Nextstrain_clade,
              values_from = clade_samples, 
              values_fill = 0)

realdata_fit <- mdl$sample(data = list(T = nrow(d), ## num observations
                                       K = ncol(d)-1,   ## num clades
                                       Y = as.matrix(d[,-1]),    ## obs
                                       N = rowSums(as.matrix(d[,-1])), ## total samples per time
                                       n_basis = 5, ## num basis funcitons
                                       B = B <- bs(1:nrow(d), df = 5, intercept=TRUE), ## design matrix of basis functions 
                                       sigma_beta = 5, ## prior sd for beta
                                       sigma_beta_rw = .5), ## prior sd for beta RW 
                           chains=4, 
                           parallel_chains=4)
## get and plot simulated pi's
draws_pi <- as_draws_df(realdata_fit$draws("pi"))

pi_est <- pivot_longer(draws_pi, cols = starts_with("pi")) %>%
  rename(chain=.chain) %>%
  group_by(name) %>%
  summarize(med = median(value),
            q10 = quantile(value, probs=0.1),
            q90 = quantile(value, probs=0.9)) %>%
  tidyr::extract(name, c("week", "clade"), "([[:alnum:]]+),([[:alnum:]]+)", remove=FALSE, convert=TRUE)

model_ests <- ggplot(pi_est, aes(x=week)) +
  geom_line(aes(y=med, color=factor(clade))) + 
  geom_ribbon(aes(fill=factor(clade), ymin=q10, ymax=q90), alpha=.3)+
  ylim(c(0,1)) +
  ylab("clade probabilities")

data_plot <- raw_d %>%
  group_by(epidate) %>%
  mutate(total_samples = sum(clade_samples)) %>%
  ungroup() %>%
  mutate(clade_prop = clade_samples / total_samples) %>%
  ggplot(aes(x=epidate, y=clade_prop, color=Nextstrain_clade)) +
  geom_line()

cowplot::plot_grid(model_ests, data_plot, nrow=2)

