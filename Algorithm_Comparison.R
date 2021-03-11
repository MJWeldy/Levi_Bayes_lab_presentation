
library(knitr)
library(reticulate) # greta uses tensorflow which is accessed through reticulate
source("install_greta.R") #only install if you are familiar with reticulate or dont use it at all greta will not run without it!!!
reticulate::use_condaenv("r-reticulate")
reticulate::py_config()
library(DiagrammeR) #Take a long time to install
library(igraph) #Takes a long time to install

library(coda)
library(bayesplot)

#MCMC engines
library(R2jags)
library(rstan)
library(nimble)
library(greta)

# Huggins Single Site --------------------------------------
# Simulating data
N <- 150 # True population size
n_occ <- 4 # Number of trapping occasions
p <- 0.50 # Probability of first detection
true_detections <- array(NA, dim = c(N, n_occ))
for (t in 1:n_occ) {
  true_detections[, t] <- rbinom(n = N, size = 1, prob = p)
}
observed <- true_detections[apply(true_detections, 1, max) == 1, ] #observed CH
(MNKA <- nrow(observed))

# JAGS --------------------------------------------------------------------

data <- list(
  y = observed,
  MNKA = MNKA,
  n_occ = n_occ
)
model_string <- textConnection(
  "
    model {
    # Likelihood
    for(i in 1:MNKA) {
      # Observation model
      for(j in 1:n_occ) {
        y[i, j] ~ dbern(p)
      } #j
    } #i
    
    # Priors
    p ~ dunif(0, 1) # Uninformative prior
    
    # Derived values
    for(t in 1:n_occ){
      p_un[t] <- (1-p)
    } #t
    N <- (MNKA / (1-prod(p_un[])))
  }"
)
parameters <- c("p", "N")
  inits <- function() {
    list()
  }
ni <- 10000 ; nt <- 1 ; nb <- 5000; nc <- 3
model1 <- jags(data, inits, parameters, model_string, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
model1_mcmc <- as.mcmc(model1)
mcmc_combo(model1_mcmc, pars = c("N", "p"))
round(model1$BUGSoutput$summary, 2)

# NIMBLE ------------------------------------------------------------------

n_data <- list(
  y = observed,
  MNKA = MNKA
)
n_constants <- list(
  n_occ = n_occ
)
Nimble_Code <- nimbleCode({
  # Likelihood
  for(i in 1:MNKA) {
    # Observation model
    for(j in 1:n_occ) {
      y[i, j] ~ dbern(p)
    } #j
  } #i
  
  # Priors
  p ~ dunif(0, 1) # Uninformative prior
  
  # Derived values
  for(t in 1:n_occ){
    p_un[t] <- (1-p)
  } #t
  N <- (MNKA / (1-prod(p_un[1:n_occ]))) #The only difference in this model is here declaring dimensions
})

n_params <- c("p", "N")
n_inits <- list( 
  
)
Nimble_Model1 <- nimbleModel(
  code = Nimble_Code,
  constants = n_constants,
  data = n_data,
  inits = n_inits
)
MCMC_Model1 <- configureMCMC(Nimble_Model1, monitors = n_params, print = T, enableWAIC = F)
Model1_MCMC <- buildMCMC(MCMC_Model1)
Comp_Model1 <- compileNimble( Nimble_Model1, showCompilerOutput = TRUE )
Comp_Model1 <- compileNimble(Model1_MCMC, project = Nimble_Model1)

niter=10000
Model1_samples <- runMCMC(Comp_Model1, niter = niter, nburnin=niter/2,nchains=3,summary = TRUE)
mcmc_combo(Model1_samples$samples, pars = c("N", "p"))
round(Model1_samples$summary$all.chains,2)

# Stan --------------------------------------------------------------------

data <- list(
  y = observed,
  MNKA = MNKA,
  n_occ = n_occ
)
stan_model <- "
  data {
    int<lower=0> MNKA;
    int<lower=0> n_occ;
    int<lower=0,upper=1> y[MNKA, n_occ];
  }
  parameters {
    real<lower=0, upper=1> p;
  }
  model {  
   for(i in 1:MNKA)
      for(j in 1:n_occ)
        target += bernoulli_lpmf(y[i, j] | p); //i and j
  }
  generated quantities {
    real pstar = (1-(1-p)^n_occ);
    real N = MNKA / pstar;
  }
  "
nc <- 4
stan1.samples <- stan(model_code = stan_model, data = data, iter = 10000, chains = nc, cores = nc, open_progress = FALSE)
mcmc_combo(as.array(stan1.samples), pars = c("N", "p"))
round(summary(stan1.samples)$summary, 2)

# greta -------------------------------------------------------------------

#capture_vec <- as_data(observed)
capture_vec <- unlist(observed)
# priors
p_greta <- beta(1, 1)
# likelihood
distribution(capture_vec) <- bernoulli(p_greta)
# derived parameters
pstar <- 1 - (1 - p_greta)^n_occ
N_hat <- MNKA / pstar

m1 <- model(p_greta, N_hat, pstar) # defining the model
plot(m1) # DAG
draws1 <- greta::mcmc(m1, n_samples = 10000, warmup=5000, n_cores=nc) # sampling

mcmc_combo(draws1, pars = c("N_hat", "p_greta"))
summary(draws1)

# Huggins p(.)-----------------------------------

# Simulating data
n_sites <- 2
N[1] <- 150 # True population size
N[2] <- 100
n_occ <- 6 # Number of trapping occasions
p <- 0.30 # Probability of first detection
true_detections <- array(NA, dim = c(N[1], n_occ, n_sites))
dim(true_detections)
for (s in 1:2) {
  for (t in 1:n_occ) {
    tmp <- rbinom(n = N[s], size = 1, prob = p)
    true_detections[1:N[s], t, s] <- tmp
  }
}
observed1 <- true_detections[apply(true_detections[,,1], 1, max) == 1, ,1]
observed2 <- true_detections[apply(true_detections[,,2], 1, max) == 1, ,2]
observed2 <- observed2[complete.cases(observed2),]

MNKA[1] <- nrow(observed1)
MNKA[2] <- nrow(observed2)
MNKA

observed <- array(NA, dim = c(MNKA[1], n_occ, n_sites))
observed[,,1] <- observed1
observed[1:MNKA[2],,2] <- observed2

# JAGS --------------------------------------------------------------------

data <- list(
  y = observed,
  n_sites = n_sites,
  MNKA = MNKA,
  n_occ = n_occ
)
model_string <- textConnection(
  "
    model {
    # Likelihood
    for(s in 1:n_sites) {
      for(i in 1:MNKA[s]) {
        # Observation model
        for(j in 1:n_occ) {
          y[i, j, s] ~ dbern(p)
        } #j
      } #i
    } #s 
    
    # Priors
    p ~ dunif(0, 1) # Uninformative prior
    
    # Derived values
    for(t in 1:n_occ){
      p_un[t] <- (1-p)
    } #t
    for (s in 1:n_sites) {
      N[s] <- (MNKA[s] / (1-prod(p_un[])))
    } #s
  }"
)
parameters <- c("p", "N")
  inits <- function() {
    list()
  }
ni <- 10000 ; nt <- 1 ; nb <- 5000; nc <- 3
model2 <- jags(data, inits, parameters, model_string, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
model2_mcmc <- as.mcmc(model2)
mcmc_combo(model2_mcmc, pars = c("N[1]","N[2]", "p"))
round(model2$BUGSoutput$summary, 2)

# NIMBLE ------------------------------------------------------------------

n_data <- list(
  y = observed,
  MNKA = MNKA
)
n_constants <- list(
  n_sites = n_sites,
  n_occ = n_occ
)
Nimble_Code <- nimbleCode({
  # Likelihood
  for(s in 1:n_sites) {
    for(i in 1:MNKA[s]) {
      # Observation model
      for(j in 1:n_occ) {
        y[i, j, s] ~ dbern(p)
      } #j
    } #i
  } #s 
  
  # Priors
  p ~ dunif(0, 1) # Uninformative prior
  
  # Derived values
  for(t in 1:n_occ){
    p_un[t] <- (1-p)
  } #t
  for (s in 1:n_sites) {
    N[s] <- (MNKA[s] / (1-prod(p_un[1:n_occ])))
  } #s
})

n_params <- c("p", "N")
n_inits <- list( 
  
)
Nimble_Model2 <- nimbleModel(
  code = Nimble_Code,
  constants = n_constants,
  data = n_data,
  inits = n_inits
)
MCMC_Model2 <- configureMCMC(Nimble_Model2, monitors = n_params, print = T, enableWAIC = F)
Model2_MCMC <- buildMCMC(MCMC_Model2)
Comp_Model2 <- compileNimble( Nimble_Model2, showCompilerOutput = TRUE )
Comp_Model2 <- compileNimble(Model2_MCMC, project = Nimble_Model2)

niter=10000
Model2_samples <- runMCMC(Comp_Model2, niter = niter, nburnin=niter/2,nchains=3,summary = TRUE)
mcmc_combo(Model2_samples$samples, pars = c("N[1]","N[2]", "p"))
round(Model2_samples$summary$all.chains,2)

# Stan --------------------------------------------------------------------

observed_stan <- observed
observed_stan[is.na(observed_stan)] <- 0
data <- list(
  y = observed_stan,
  n_sites = n_sites,
  M = MNKA[1],
  MNKA = MNKA,
  n_occ = n_occ
)

stan_model <- "
  data {
    int<lower=0> n_sites;
    int<lower=0> n_occ;
    int<lower=0> M;
    int<lower=0> MNKA[n_sites];
    int<lower=0,upper=1> y[M, n_occ, n_sites];
  }
  parameters {
    real<lower=0, upper=1> p;
  }
  model {  
   for(s in 1:n_sites)
    for(i in 1:MNKA[s])
      for(j in 1:n_occ)
        target += bernoulli_lpmf(y[i, j, s] | p); //s i j
  }
  generated quantities {
    real N[n_sites];
    real<lower=0> pstar;
    
    pstar = (1-(1-p)^n_occ);
    for(s in 1:n_sites)
      N[s] = MNKA[s] / pstar; // s
  }
  "
nc <- 4
stan2.samples <- stan(model_code = stan_model, data = data, iter = 10000, chains = nc, warmup=5000, cores = nc, open_progress = FALSE)
mcmc_combo(as.array(stan2.samples), pars = c("N[1]","N[2]", "p"))
round(summary(stan2.samples)$summary, 2)

# greta -------------------------------------------------------------------

capture_vec <- unlist(unlist(rbind(observed1,observed2)))
capture_vec <- as.vector(t(capture_vec))
p_det <- as.vector(c(rep(1,length(observed1)),rep(2,length(observed2))))

# priors
p_greta <- beta(1, 1, dim = 2)

# likelihood
distribution(capture_vec) <- bernoulli(p_greta[p_det])
# derived parameters
pstar <- 1 - (1 - p_greta)^6
N_hat <- MNKA / pstar

m2 <- model(p_greta, N_hat, pstar) # defining the model
plot(m2) # DAG
draws2 <- greta::mcmc(m2, n_samples = 10000, warmup=5000, n_cores=nc) # sampling
mcmc_combo(draws2, pars = c("N_hat[1,1]","N_hat[2,1]","p_greta[1,1]",  "p_greta[2,1]"))
summary(draws2)

# Huggins p(site)-----------------------------------
# Simulating data
n_sites <- 2
N[1] <- 150 # True population size
N[2] <- 100
n_occ <- 6 # Number of trapping occasions
p[1] <- 0.30 # Probability of first detection
p[2] <- 0.50 # Probability of first detection
true_detections <- array(NA, dim = c(N[1], n_occ, n_sites))
dim(true_detections)
for (s in 1:2) {
  for (t in 1:n_occ) {
    tmp <- rbinom(n = N[s], size = 1, prob = p[s])
    true_detections[1:N[s], t, s] <- tmp
  }
}
observed1 <- true_detections[apply(true_detections[,,1], 1, max) == 1, ,1]
observed2 <- true_detections[apply(true_detections[,,2], 1, max) == 1, ,2]
observed2 <- observed2[complete.cases(observed2),]

MNKA[1] <- nrow(observed1)
MNKA[2] <- nrow(observed2)
MNKA

observed <- array(NA, dim = c(MNKA[1], n_occ, n_sites))
observed[,,1] <- observed1
observed[1:MNKA[2],,2] <- observed2

# JAGS --------------------------------------------------------------------

data <- list(
  y = observed,
  n_sites = n_sites,
  MNKA = MNKA,
  n_occ = n_occ
)
model_string <- textConnection(
  "
    model {
    # Likelihood
    for(s in 1:n_sites) {
      for(i in 1:MNKA[s]) {
        # Observation model
        for(j in 1:n_occ) {
          y[i, j, s] ~ dbern(p[s])
        } #j
      } #i
    } #s 
    
    # Priors
    for(s in 1:n_sites) {
      p[s] ~ dunif(0, 1) # Uninformative prior
    } #s
    
    # Derived values
    for (s in 1:n_sites) {
      for(t in 1:n_occ){
        p_un[s,t] <- (1-p[s])
      } #t
      N[s] <- (MNKA[s] / (1-prod(p_un[s,])))
    } #s
  }"
)
parameters <- c("p", "N")
inits <- function() {
  list()
}
ni <- 10000 ; nt <- 1 ; nb <- 5000; nc <- 3
model3 <- jags(data, inits, parameters, model_string, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
model3_mcmc <- as.mcmc(model3)
mcmc_combo(model3_mcmc, pars = c("N[1]","N[2]", "p[1]","p[2]"))
round(model3$BUGSoutput$summary, 2)

# Stan --------------------------------------------------------------------
observed_stan <- observed
observed_stan[is.na(observed_stan)] <- 0
data <- list(
  y = observed_stan,
  n_sites = n_sites,
  M = MNKA[1],
  MNKA = MNKA,
  n_occ = n_occ
)

stan_model <- "
  data {
    int<lower=0> n_sites;
    int<lower=0> n_occ;
    int<lower=0> M;
    int<lower=0> MNKA[n_sites];
    int<lower=0,upper=1> y[M, n_occ, n_sites];
  }
  parameters {
    real<lower=0, upper=1> p[n_sites];
  }
  model {  
   for(s in 1:n_sites)
    for(i in 1:MNKA[s])
      for(j in 1:n_occ)
        target += bernoulli_lpmf(y[i, j, s] | p[s]); //s i j
  }
  generated quantities {
    real N[n_sites];
    real<lower=0> pstar[n_sites];
    
    for(s in 1:n_sites) {
      pstar[s] = (1-(1-p[s])^n_occ);
      N[s] = MNKA[s] / pstar[s];
    } //s
  }
  "
nc <- 4
stan3.samples <- stan(model_code = stan_model, data = data, iter = 10000, chains = nc, warmup=5000, cores = nc, open_progress = FALSE)
mcmc_combo(as.array(stan3.samples), pars = c("N[1]","N[2]", "p[1]", "p[2]"))
round(summary(stan3.samples)$summary, 2)


# greta -------------------------------------------------------------------

capture_vec <- unlist(unlist(rbind(observed1,observed2)))
capture_vec <- as.vector(t(capture_vec))
p_det <- as.vector(c(rep(1,length(observed1)),rep(2,length(observed2))))

# priors
p_greta <- beta(1, 1, dim = 2)
# likelihood
distribution(capture_vec) <- bernoulli(p_greta[p_det])
# derived parameters
pstar <- 1 - (1 - p_greta)^6
N_hat <- MNKA / pstar

m3 <- model(p_greta, N_hat, pstar) # defining the model
plot(m3) # DAG
draws3 <- greta::mcmc(m3, n_samples = 10000, warmup=5000, n_cores=nc) # sampling
mcmc_combo(draws3, pars = c("N_hat[1,1]","N_hat[2,1]","p_greta[1,1]",  "p_greta[2,1]"))
summary(draws3)

# CJS ---------------------------------------------------------------------
# Simulating data
n_occ <- 6                   # Number of capture occasions
marked <- rep(50, n_occ-1)   # Annual number of newly marked individuals
phi <- rep(0.65, n_occ-1)
p <- rep(0.4, n_occ-1)
# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n_occ-1, nrow = sum(marked))
P <- matrix(p, ncol = n_occ-1, nrow = sum(marked))
CH <- matrix(0, ncol = n_occ, nrow = sum(marked))
marking_occ <- rep(1:length(marked), marked[1:length(marked)])
# Fill the CH matrix
i<-1
for (i in 1:sum(marked)){
  CH[i, marking_occ[i]] <- 1       # Write an 1 at the release occasion
  if (marking_occ[i]==n_occ) next
  for (t in (marking_occ[i]+1):n_occ){
    survive_occasion <- rbinom(1, 1, PHI[i,t-1])
    if (survive_occasion==0) break 
    rp <- rbinom(1, 1, P[i,t-1])
    if (rp==1) CH[i,t] <- 1
  } #t
} #i
get.first.capture <- function(x) min(which(x!=0))
first_capture <- apply(CH, 1, get.first.capture)
get.last.capture <- function(x) max(which(x!=0))
last_capture <- apply(CH, 1, get.last.capture)

# JAGS --------------------------------------------------------------------
data <- list(
  y = CH,
  n_ind = dim(CH)[1],
  n_occ = dim(CH)[2],
  f = first_capture
)
model_string <- textConnection(
  "
  model {
  
  # Likelihood
  for (i in 1:n_ind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n_occ){
      # State process
      z[i,t] ~ dbern(phi * z[i,t-1])
      # Observation process
      y[i,t] ~ dbern(mu[i,t])
      mu[i,t] <- p * z[i,t]
    } #t
  } #i
  
  # Priors and Constraints
  phi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  # Derived values
  }
"
)
parameters <- c("p","phi")
cjs.z.init <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}
set_initial_value <- function() {
  list( 
    z = cjs.z.init(CH)  
  )
}
ni <- 10000 ; nt <- 1 ; nb <- 5000 ; nc <- 3
model4 <- jags(data, set_initial_value, parameters, model_string, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
model4_mcmc <- as.mcmc(model4)
mcmc_combo(model4_mcmc, pars = c("phi", "p"))
round(model4$BUGSoutput$summary, 2)

# Stan -------------------------------------------------------------------
data <- list(
  CH = CH,
  n_ind = dim(CH)[1],
  n_occ = dim(CH)[2],
  f = first_capture,
  l = last_capture
)
stan_model <- "
  /**
   * following section 1.2.1 of:
   * http://www.maths.otago.ac.nz/home/resources/theses/PhD_Matthew_Schofield.pdf
   * https://discourse.mc-stan.org/t/cjs-log-likelihood/15112
   */
  data {
    int<lower=2> n_occ;                      // capture events
    int<lower=0> n_ind;                      // number of individuals
    int<lower=0, upper=n_occ+1> f[n_ind];     // f[i]: ind i first capture
    int<lower=0, upper=n_occ+1> l[n_ind];     // l[i]:  ind i last capture
    int<lower=0,upper=1> CH[n_ind,n_occ];    // CH[i,k]: individual i captured at k
  }
  transformed data {
    int<lower=0,upper=n_ind> n_captured[n_occ];  // n_capt[k]: num captured at k
    n_captured = rep_array(0,n_occ);
    for (i in 1:n_ind)
      for (k in 1:n_occ)
        n_captured[k] = n_captured[k] + CH[i,k]; //i k 
  }
  parameters {
    real<lower=0,upper=1> phi;  // phi[k]: Pr[alive at k + 1 | alive at k]
    real<lower=0,upper=1> p;      // p[k]: Pr[capture at k]
  }
  transformed parameters {
    vector<lower=0,upper=1>[n_occ] chi;   // chi[k]: Pr[no capture >  k | alive at k]
    vector[n_ind] log_lik;
    {
      int k;
      chi[n_occ] = 1.0;              
      k = n_occ - 1;
      while (k > 0) {
        chi[k] = (1 - phi) + phi * (1 - p) * chi[k+1];
        k = k - 1;
      }
    }
    
    for (i in 1:n_ind) {
      log_lik[i] = 0;
      if (l[i] > 0) {
        for (k in (f[i]+1):l[i]) {
          log_lik[i] +=log(phi);     // i survived from k-1 to k
          if (CH[i,k] == 1)
            log_lik[i] +=log(p);       // i captured at k
          else
            log_lik[i] +=log1m(p);     // i not captured at k
        }
        log_lik[i] +=log(chi[l[i]]);   // i not seen after last[i]
      }
    }
  }
  model {
    target += sum(log_lik);
  }
  generated quantities {
    // phi[K-1] and p(K) not identified, but product is
    real beta;
    vector<lower=0>[n_occ] pop_hat;  // population
    
    beta = phi * p;
    for (k in 1:n_occ)
      pop_hat[k] = n_captured[k] / p;  // k
  }
"
inits <- function() list(phi = runif(1, 0, 1),
                         p = runif(1, 0, 1))
## MCMC settings
nc <- 4
stan4.samples <- stan(model_code = stan_model, data = data, iter = 10000, chains = nc, warmup=5000, cores = nc, open_progress = FALSE)
mcmc_combo(as.array(stan4.samples), pars = c("phi", "p"))
round(summary(stan4.samples)$summary, 2)

# greta -------------------------------------------------------------------
obs_id <- apply(CH, 1, function(x) seq(min(which(x > 0)), max(which(x > 0)), by = 1)[-1])
obs_id <- unlist(obs_id)
capture_vec <- apply(CH, 1, function(x) x[min(which(x > 0)):max(which(x > 0))][-1])
capture_vec <- unlist(capture_vec)
# dummy variables
alive_data <- ones(length(obs_id))            # definitely alive
not_seen_last <- last_capture != n_occ              # ignore observations in last timestep
final_observation <- ones(sum(not_seen_last)) # final observation

#capture_vec <- as_data(observed)

# priors
phi <- beta(1, 1, dim = 1)
p <- beta(1, 1, dim = 1)
# derived parameter
chi <- ones(n_occ)
for (i in seq_len(n_occ - 1)) {
  tn <- n_occ - i
  chi[tn] <- (1 - phi) + phi * (1 - p) * chi[tn + 1]
}
# likelihood
distribution(alive_data) <- bernoulli(phi)
distribution(capture_vec) <- bernoulli(p)
distribution(final_observation) <- bernoulli(chi[last_capture[not_seen_last]])
# defining the model
m4 <- model(phi,p) #objects to sample
# sampling
draws4 <- greta::mcmc(m4, n_samples = 10000, warmup=5000, n_cores=nc)
plot(m4)
mcmc_combo(draws4, pars = c("phi","p"))
summary(draws4)
