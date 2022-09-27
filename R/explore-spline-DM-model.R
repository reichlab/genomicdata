## simulate clade data from Stan spline model

library(tidyverse)
library(cmdstanr)
library(posterior)
library(splines)

theme_set(theme_bw())

mdl <- cmdstan_model(stan_file = "R/simple-dirichlet-4.stan")

########
## simulate fake data
T <- 100
K <- 3
t <- 1:T
n_basis <- max(round(T/10), 4)
B <- bs(t, df = n_basis, intercept=TRUE)
N <- rep(200, times=T)  

# simulate betas
mu_beta <- 0
sigma_beta <- 5
sigma_beta_rw <- 0.2
sim_betas_cauchy <- function(n_basis, K, s_beta, s_beta_rw) {
  beta <- matrix(NA, nrow=n_basis, ncol=K)
  beta[1,] <- rnorm(K, 0, s_beta)
  for(j in 2:n_basis) {
    beta[j,] <- beta[j-1, ] + rcauchy(K, 0, s_beta_rw)
  }
  return(beta)
}
sim_betas_norm <- function(n_basis, K, mu_beta, s_beta) {
  beta <- matrix(rnorm(n_basis*K, mu_beta, s_beta), nrow=n_basis, ncol=K)
  return(beta)
}
sim_betas_halfnormal <- function(n_basis, K, mu_beta, s_beta) {
  beta <- matrix(abs(rnorm(n_basis*K, mu_beta, s_beta)), nrow=n_basis, ncol=K)
  return(beta)
}
sim_betas_halfcauchy <- function(n_basis, K, mu_beta, s_beta) {
  beta <- matrix(abs(rcauchy(n_basis*K, mu_beta, s_beta)), nrow=n_basis, ncol=K)
  return(beta)
}


nsim <- 20
betas <- replicate(nsim, 
                   sim_betas_halfcauchy(n_basis=n_basis, K=K, mu_beta = mu_beta, s_beta=sigma_beta), 
                   simplify=FALSE)

phi_plots <- vector(mode="list", length=nsim)
prob_plots <- vector(mode="list", length=nsim)
#par(mfrow=c(4, 5), mar=c(.2, .2, .2, .2))
y <- array(NA, dim=c(T, K, nsim))
for(i in 1:nsim){
  # simulate ys
  phi <- B %*% betas[[i]]
  for(t in 1:T)
    y[t,,i] <- t(rmultinom(1, 
                           size = N[1], 
                           prob = gtools::rdirichlet(1, phi[t,])))
  
  #plot(1:T, y[,1,i], type="l", ylim=c(0,N[1]), xaxt="n", yaxt="n")
  #points(1:T, y[,2,i], type="l", col="red")
  #points(1:T, y[,3,i], type="l", col="blue")


  ## run HMC for one beta and compare to truth
  
  mdl_fit <- mdl$sample(data = list(T = T, ## num observations
                                    K = K,   ## num clades
                                    Y = y[,,i],    ## obs
                                    n_basis = n_basis, ## num basis funcitons
                                    B = B,        ## design matrix of basis functions 
                                    sigma_beta = 5, ## prior sd for beta
                                    mu_beta = 0), ## prior mean for beta
                        chains=4, 
                        parallel_chains=4)
  
  ## get and plot simulated pi's
  draws_phi <- as_draws_df(mdl_fit$draws("phi"))
  
  phi_est <- pivot_longer(draws_phi, cols = starts_with("phi")) %>%
    rename(chain=.chain) %>%
    tidyr::extract(name, c("week", "clade"), "([[:alnum:]]+),([[:alnum:]]+)", remove=FALSE, convert=TRUE) %>% 
    group_by(week, chain, .iteration, .draw) %>%
    mutate(total_phi = sum(value)) %>%
    ungroup() %>%
    mutate(exp_prob = value/total_phi) %>%
    group_by(name) %>%
    summarize(med = median(value),
              q10 = quantile(value, probs=0.1),
              q90 = quantile(value, probs=0.9),
              med_prob = median(exp_prob),
              q10_prob = quantile(exp_prob, probs=0.1),
              q90_prob = quantile(exp_prob, probs=0.9))  %>%
    tidyr::extract(name, c("week", "clade"), "([[:alnum:]]+),([[:alnum:]]+)", remove=FALSE, convert=TRUE)
  
  true_phis <- B %*% betas[[i]] |> 
    as_tibble() |> 
    rename(phi_1 = V1,
           phi_2 = V2,
           phi_3 = V3) |> 
    select(starts_with("phi")) |> 
    mutate(week = 1:n()) %>%
    pivot_longer(cols = starts_with("phi"), names_to = "clade", values_to = "phi") %>%
    group_by(week) %>%
    mutate(exp_prob = phi/sum(phi))
  
  dat <- data.frame(y[,,i]) %>%
    mutate(week=1:n()) %>%
    pivot_longer(cols=-week) %>%
    mutate(clade = substr(name, 2, 3)) %>%
    group_by(week) %>%
    mutate(freq = value/sum(value))
  
  ## plot phi
  phi_plots[[i]] <- ggplot(phi_est, aes(x=week)) +
    geom_line(aes(y=med, color=factor(clade))) + 
    geom_line(data=true_phis, aes(y=phi, group=clade)) +
    geom_ribbon(aes(fill=factor(clade), ymin=q10, ymax=q90), alpha=.3)+
    ylab("clade dirichlet params") +
    theme(legend.position = "none")
  
  ## plot probs
  prob_plots[[i]] <- ggplot(phi_est, aes(x=week)) +
    geom_line(aes(y=med_prob, color=factor(clade))) + 
    geom_line(data=true_phis, aes(y=exp_prob, group=clade)) +
    geom_point(data=dat, aes(y=freq, color=clade), alpha=.2) +
    geom_ribbon(aes(fill=factor(clade), ymin=q10_prob, ymax=q90_prob), alpha=.3)+
    scale_y_continuous("clade probabilities", limits=c(0,1)) +
    theme(legend.position = "none")
}  

phips <- gridExtra::marrangeGrob(phi_plots, nrow=4, ncol=5)
ggsave("phi-plots.pdf", phips, width = 12, height=8)

prps <- gridExtra::marrangeGrob(prob_plots, nrow=4, ncol=5)
ggsave("prob-plots.pdf", prps, width = 12, height=8)

########
## estimate model using real data

## selecting clade with >10K samples total
## and reserving the last ~year of data as a test set
clades <- c("20A", "20C", 
            "20G", "20I (Alpha, V1)", 
            "21I (Delta)", "21J (Delta)", 
            "21K (Omicron)")
raw_d <- read_csv("data/genomic_data_for_modeling.csv") %>%
  filter(Nextstrain_clade %in% clades, 
         epidate <= as.Date("2021-06-01"),
         epidate > as.Date("2020-04-01")) 

d <- raw_d %>%
  pivot_wider(id_cols = epiweek_year, 
              names_from = Nextstrain_clade,
              values_from = clade_samples, 
              values_fill = 0)

realdata_fit <- mdl$sample(data = list(T = nrow(d), ## num observations
                                       K = ncol(d)-1,   ## num clades
                                       Y = as.matrix(d[,-1]),    ## obs
                                       n_basis = 8, ## num basis funcitons
                                       B = B <- bs(1:nrow(d), df = 8, intercept=TRUE), ## design matrix of basis functions 
                                       mu_beta = 0, ## prior mean for beta
                                       sigma_beta = 5), ## prior sd for beta RW 
                           chains=4, 
                           parallel_chains=4)
# get and plot simulated pi's
draws_phi <- as_draws_df(realdata_fit$draws("phi"))

phi_est <- pivot_longer(draws_phi, cols = starts_with("phi")) %>%
  rename(chain=.chain) %>%
  tidyr::extract(name, c("week", "clade"), "([[:alnum:]]+),([[:alnum:]]+)", remove=FALSE, convert=TRUE) %>% 
  group_by(week, chain, .iteration, .draw) %>%
  mutate(total_phi = sum(value)) %>%
  ungroup() %>%
  mutate(exp_prob = value/total_phi) %>%
  group_by(name) %>%
  summarize(med = median(value),
            q10 = quantile(value, probs=0.1),
            q90 = quantile(value, probs=0.9),
            med_prob = median(exp_prob),
            q10_prob = quantile(exp_prob, probs=0.1),
            q90_prob = quantile(exp_prob, probs=0.9))  %>%
  tidyr::extract(name, c("week", "clade"), "([[:alnum:]]+),([[:alnum:]]+)", remove=FALSE, convert=TRUE) %>%
  mutate(clade = factor(clade, levels = 1:7, labels=colnames(d)[-1]))

data_to_plot <- raw_d %>%
  filter(Nextstrain_clade %in% clades) %>%
  group_by(epidate) %>%
  mutate(total_samples = sum(clade_samples)) %>%
  ungroup() %>%
  mutate(clade_prop = clade_samples / total_samples) %>%
  mutate(clade = factor(Nextstrain_clade, levels=levels(phi_est$clade))) %>%
  group_by(clade) %>%
  arrange(epidate) %>%
  mutate(week = as.numeric(difftime(epidate, as.Date("2020-03-29"), units = "week")))
  

ggplot(phi_est, aes(x=week)) +
  geom_point(data=data_to_plot, aes(y=clade_prop, color=factor(clade)), alpha=.3) +
  geom_line(aes(y=med_prob, color=factor(clade))) + 
  geom_ribbon(aes(fill=factor(clade), ymin=q10_prob, ymax=q90_prob), alpha=.3) +
  ylim(c(0,1)) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  ylab("clade probabilities")




