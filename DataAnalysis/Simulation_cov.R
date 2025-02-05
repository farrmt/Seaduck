library(nimble)
library(coda)
library(parallel)
library(tidyverse)

#Logit function
logit <- function(pp) 
{ 
  log(pp) - log(1-pp) 
}

#Inverse logit
expit <- function(eta) 
{
  1/(1+exp(-eta))
}

#Generate miss ID rates
missID.fun <- function(nspecies, nmissID)
{
  phi.pi <- matrix(0, ncol = nspecies, nrow = nspecies)
  phi.pi[sample(x = (1:nspecies^2)[-seq(1, 100, nspecies+1)], size = nmissID)] <- runif(nmissID, 0, 0.5/nspecies)
  diag(phi.pi) <- 1 - apply(phi.pi, MARGIN = 1, sum)
  # diag(phi.pi) <- runif(nspecies, 0.8, 1)
  # phi.pi <- phi.pi/apply(phi.pi, MARGIN = 1, sum)
  return(phi.pi)
}

#Simulate values between 1 - 2 SD from mean
cov.cam.fun <- function(nreps = 10){
  valrange <- data.frame(x1 = seq(1,2 - 0.001,0.001), 
                         x2 = seq(1 + 0.001, 2, 0.001), 
                         value = seq(1.0005, by = 0.001, length.out = 1000))
  valrange <- valrange %>% mutate(prob = (pnorm(x2) - pnorm(x1))/(pnorm(2) - pnorm(1)))
  value <- sample(valrange$value, size = nreps, replace = T, prob = valrange$prob)
  return(value)
}


#Number of years
nyears <- 10

#Number of sites
nsites <- 75

#Number of species
nspecies <- 10

#Simulation iteration start
#for(iter in 21:100){

#Expected population growth rate intercept
mu.g0 <- log(runif(1, 0.95, 1.05))
sd.g0 <- runif(1, 0.1, 0.25)
gamma0 <- rnorm(nspecies, mu.g0, sd.g0)
mu.g0 <- mean(gamma0)
sd.g0 <- sd(gamma0)

#Expected population growth rate covariate effect
# mu.g1 <-  runif(1, -0.1, 0.1)
mu.g1 <-  runif(1, 0.25, 0.5)
# sd.g1 <- runif(1, 0.1, 0.5)
sd.g1 <- runif(1, 0.1, 0.25)
gamma1 <- rnorm(nspecies, mu.g1, sd.g1)
mu.g1 <- mean(gamma1)
sd.g1 <- sd(gamma1)

#Expected population growth rate covariate
cov.t <- rnorm(nyears - 1, mean = 0, sd = 1)

#Expected population growth rate
gamma <- matrix(data = NA, nrow = nspecies, ncol = nyears - 1)

for(i in 1:nspecies){
  gamma[i,] <- exp(gamma0[i] + gamma1[i] * cov.t)
}

#Expected population abundance intercept in year 1 at the site-level
mu.b0 <- runif(1, 1, 2)
# sd.b0 <- runif(1, 1, 2)
# mu.b0 <- runif(1, 0.5, 1)
sd.b0 <- runif(1, 1, 2)
beta0 <- rnorm(nspecies, mu.b0, sd.b0)
mu.b0 <- mean(beta0)
sd.b0 <- sd(beta0)

#Spatiotemporal covariate effects
# mu.b1 <- runif(1, -0.5, 0.5)
# mu.b1 <- runif(1, 0.25, 0.5)
mu.b1 <- runif(1, 0.5, 1)
# sd.b1 <- runif(1, 1, 2)
sd.b1 <- runif(1, 0.1, 0.25)
beta1 <- rnorm(nspecies, mu.b1, sd.b1)
mu.b1 <- mean(beta1)
sd.b1 <- sd(beta1)

#Spatiotemporal covariate
cov.st <- matrix(NA, nrow = nyears, ncol = nsites)

for(t in 1:nyears){
  cov.st[t,1:nsites] <- rnorm(nsites, mean = 0, sd = 1)
  cov.st[t,1:nsites] <- scale(cov.st[t,1:nsites])
}

#Expected population abundance at the site-level
lambda <- array(NA, dim = c(nspecies, nyears, nsites))

#Latent site-level population abundance
n <- array(NA, dim = c(nspecies, nyears, nsites))

#Generate values for the above
for(i in 1:nspecies){
  lambda[i,1,1:nsites] <- exp(beta0[i] + beta1[i] * cov.st[1,1:nsites])
  n[i,1,1:nsites] <- rpois(nsites, lambda[i,1,1:nsites])
  for(t in 2:nyears){
    lambda[i,t,1:nsites] <- exp(beta0[i] + beta1[i] * cov.st[t,1:nsites] + log(prod(gamma[i,1:(t-1)])))
    n[i,t,1:nsites] <- rpois(nsites, lambda[i,t,1:nsites])
  }
}

#-Observation process-#

#Detection probability
p <- runif(1, 0.1, 1)

#Movement rate
alpha <- runif(1, 0.5, 1)

#Expected miss ID rate
phi.pi <- missID.fun(nspecies, nmissID = 15)

#Observed counts
y <- array(NA, dim = c(nspecies, nyears, nsites))

for(t in 1:nyears){
  for(j in 1:nsites){
    tmp <- NULL
    for(i in 1:nspecies){
      tmp <- cbind(tmp, rmultinom(1, n[i,t,j], phi.pi[i,]))
    }#end i
    y[,t,j] <- rbinom(n = nspecies, apply(tmp, MARGIN = 1, sum), prob = p * alpha)
  }#end j
}#end t

#-Camera data-#

#Number of camera replicates
nreps <- 10

# #Covariate data
cov.cam <- cov.cam.fun(nreps = nreps)

#Species-specific abundance
# theta <- exp(beta0)
theta <- array(NA, dim = c(nspecies, nreps))

for(k in 1:nreps){
  for(i in 1:nspecies){
    theta[i,k] <- exp(beta0[i] + beta1[i] * cov.cam[k])
  }#end i
}#end k

#Community expected abundance
THETA <- apply(theta, 2, sum)
# THETA <- sum(theta)

#Community composition
pi <- exp(beta0)/sum(exp(beta0))
# pi <- theta/THETA

#Front facing and point of view camera data
FF <- array(NA, dim = c(nspecies, nreps))

#Observer data
OBS <- array(NA, dim = c(nspecies, nreps, 2))

for(k in 1:nreps){
  FF[,k] <- rpois(nspecies, theta[,k])
  # FF[,k] <- rpois(nspecies, theta) 
  tmp <- NULL 
  for(i in 1:nspecies){
    tmp <- cbind(tmp, rmultinom(1, FF[i,k], phi.pi[i,]))
  }#end i
  OBS[,k,1] <- rbinom(n = nspecies, apply(tmp, MARGIN = 1, sum), prob = p * alpha)
  OBS[,k,2] <- rbinom(n = nspecies, apply(tmp, MARGIN = 1, sum), prob = p * alpha)
}#end k

#-Extra parameters-#
epsilon.hat <- p*alpha
phi <- pi %*% phi.pi
correction <- as.numeric(epsilon.hat * phi/pi)
log.correction <- log(correction)
min.spec <- which.min(apply(apply(y,c(1,2), sum), 1, sum))
max.spec <- which.max(apply(apply(y,c(1,2), sum), 1, sum))
target <- c(min.spec, max.spec)

#-Output-#
output <- data.frame(matrix(NA, nrow = 1, ncol = 7))
colnames(output) <- c("model", "parameter.type", "parameter", "truth", "mean", "sd", "rhat")

#----------------------#
#-Auxiliary Model Code-#
#----------------------#

code <- nimbleCode({
  
  #PRIORS
  
  #Population density
  mu.b0 ~ dnorm(0, 0.01)
  tau.b0 ~ dgamma(1, 0.001)
  sd.b0 <- 1/sqrt(tau.b0)
  
  mu.b1 ~ dnorm(0, 0.01)
  tau.b1 ~ dgamma(1, 0.001)
  sd.b1 <- 1/sqrt(tau.b1)
  
  for(i in 1:nspecies){
    
    #Population density
    beta0[i] ~ dnorm(mu.b0, tau.b0) 
    beta1[i] ~ dnorm(mu.b1, tau.b1)
    
    #Correction factor
    log.correction[i] <- log(correction[i])
    correction[i] <- epsilon.hat * phi[i]/pi[i]
    phi.ones[i] <- 1
    
    #LIKELIHOOD
    
    exp.beta[i] <- exp(beta0[i])
    pi[i] <- exp(beta0[i])/sum(exp.beta[1:nspecies])
    
    for(k in 1:nreps){
      
      log(theta[i,k]) <- beta0[i] + beta1[i] * cov.cam[k]
      pi.star[i,k] <- theta[i,k]/THETA[k]
      phi.star[i,k] <- (phi[i] * sum(exp.beta[1:nspecies]) * exp(beta1[i] * cov.cam[k]))/THETA[k] 
      
    }#end k
    
  }#end i
  
  for(k in 1:nreps){
    
    FF[1:nspecies,k] ~ dmulti(pi.star[1:nspecies,k], FF.total[k])
    
    #Front facing camera total abundance
    # FF.total[k] ~ dpois(THETA)
    FF.total[k] ~ dpois(THETA[k])
    
    THETA[k] <- sum(theta[1:nspecies,k])
    
    for(o in 1:nobs){
      
      OBS[1:nspecies,k,o] ~ dmulti(phi.star[1:nspecies,k], OBS.total[k,o])
      
      # OBS.total[k,o] ~ dpois(THETA * epsilon.hat)
      OBS.total[k,o] ~ dpois(THETA[k] * epsilon.hat)
      
    }#end o
    
  }#end k
  
  #Composition of latent abundance (corrected for imperfect detection)
  phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  
  #Derived product of movement and detection
  epsilon.hat ~ dnorm(0, 0.01)
  
  #Community-wide expected abundance
  # THETA <- sum(theta[1:nspecies])
  
})

#-Compile data-#

data <- list(FF = FF,
             FF.total = apply(FF, 2, sum),
             OBS = OBS,
             OBS.total = apply(OBS, c(2,3), sum),
             cov.cam = cov.cam)

constants <- list(nspecies = nspecies, nreps = nreps, nobs = 2)

#-Initial values-#

inits <- function(){list(beta0 = beta0 + runif(nspecies, -0.1, 0.1), 
                         mu.b0 = mu.b0 + runif(1, -0.1, 0.1), 
                         sd.b0 = sd.b0 + runif(1, 0, 0.25),
                         beta1 = beta1 + runif(nspecies, -0.1, 0.1),
                         mu.b1 = mu.b1 + runif(1, -0.1, 0.1),
                         sd.b1 = sd.b1 + runif(1, 0, 0.25),
                         epsilon.hat = alpha*p + runif(1, -0.1, 0.1),
                         correction = correction + runif(1, -0.1, 0.1)
)}

#-Parameters to save-#

params <- c(
  "mu.b0",
  "sd.b0",
  "beta0",
  "mu.b1",
  "sd.b1",
  "beta1",
  "epsilon.hat",
  "pi",
  "phi"
)

params2 <- "log.correction"

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params, monitors2 = params2)

MCMCconf$removeSampler(target = c("beta0", "beta1"))

for(i in 1:nspecies){
  MCMCconf$addSampler(target = c(paste0("beta0[", i, "]"), paste0("beta1[", i, "]")),
                      type = "AF_slice")
}

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

aux.out <- runMCMC(compiled.model$MCMC,
                   niter = ni, nburnin = nb,
                   nchains = nc, thin = nt, thin2 = nt,
                   samplesAsCodaMCMC = TRUE)

#-Output-#
for(k in 1:length(params)){
  if(length(get(params[k])) == 1){
    output <- rbind(output,
                    data.frame(model = "auxiliary",
                               parameter.type = "community",
                               parameter = params[k],
                               truth = get(params[k]),
                               mean = summary(aux.out$samples)[[1]][params[k],"Mean"],
                               sd = summary(aux.out$samples)[[1]][params[k], "SD"],
                               rhat = as.numeric(coda::gelman.diag(aux.out$samples[1:nc][,params[k]])[[1]][,1])))
  }else{
    for(j in 1:2){
      output <- rbind(output,
                      data.frame(model = "auxiliary",
                                 parameter.type = c("min.spec", "max.spec")[j],
                                 parameter = params[k],
                                 truth = get(params[k])[target[j]],
                                 mean = summary(aux.out$samples)[[1]][paste0(params[k], "[", target[j], "]"),"Mean"],
                                 sd = summary(aux.out$samples)[[1]][paste0(params[k], "[", target[j], "]"),"SD"],
                                 rhat = as.numeric(coda::gelman.diag(aux.out$samples[1:nc][,paste0(params[k], "[", target[j], "]")])[[1]][,1])))
    }
  }
}

for(j in 1:2){
  output <- rbind(output,
                  data.frame(model = "auxiliary",
                             parameter.type = c("min.spec", "max.spec")[j],
                             parameter = params2,
                             truth = get(params2)[target[j]],
                             mean = summary(aux.out$samples2)[[1]][paste0(params2, "[", target[j], "]"),"Mean"],
                             sd = summary(aux.out$samples2)[[1]][paste0(params2, "[", target[j], "]"),"SD"],
                             rhat = as.numeric(coda::gelman.diag(aux.out$samples2[1:nc][,paste0(params2, "[", target[j], "]")])[[1]][,1])))
}

output <- output[-1,]

#---------------------------------------#
#-Single Data Single Species Model Code-#
#---------------------------------------#

code <- nimbleCode({
  
  #PRIORS
  
  #Population growth
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)  
  
  #Population density
  beta0 ~ dnorm(0, 0.01) 
  beta1 ~ dnorm(0, 0.01)
  
  for(j in 1:nsites){
    
    y[1,j] ~ dpois(lambda[1,j])
    
    lambda[1,j] <- exp(beta0 + beta1 * cov.st[1,j])
    
  }#end j
  
  for(t in 2:nyears){
    
    log(gamma[t-1]) <- gamma0 + gamma1 * cov.t[t-1]
    
    for(j in 1:nsites){
      
      y[t,j] ~ dpois(lambda[t,j])
      
      lambda[t,j] <- exp(beta0 + log(prod(gamma[1:(t-1)])) + beta1 * cov.st[t,j])
      
    }#end j
    
  }#end t
  
})

#-Compile data-#

data <- list(y = y[min.spec,,],
             cov.t = cov.t,
             cov.st = cov.st)

constants <- list(nyears = nyears, nsites = nsites)

#-Initial values-#

inits <- function(){list(beta0 = beta0[min.spec] + runif(1, -0.1, 0.1), 
                         beta1 = beta1[min.spec] + runif(1, -0.1, 0.1), 
                         gamma0 = gamma0[min.spec] + runif(1, -0.1, 0.1), 
                         gamma1 = gamma1[min.spec] + runif(1, -0.1, 0.1)
)}

#-Parameters to save-#

params <- c(
  "beta0",
  "beta1",
  "gamma0",
  "gamma1"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1"))

MCMCconf$addSampler(target = c("gamma0", "gamma1"),
                    type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  output <- rbind(output,
                  data.frame(model = "sdm.ss",
                             parameter.type = "min.spec",
                             parameter = params[k],
                             truth = get(params[k])[target[1]],
                             mean = summary(out)[[1]][params[k],"Mean"],
                             sd = summary(out)[[1]][params[k], "SD"],
                             rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
}

#-Compile data-#

max.spec <- which.max(apply(apply(y,c(1,2), sum), 1, sum))

data <- list(y = y[max.spec,,],
             cov.t = cov.t,
             cov.st = cov.st)

#-Initial values-#

inits <- function(){list(beta0 = beta0[max.spec] + runif(1, -0.1, 0.1), 
                         beta1 = beta1[max.spec] + runif(1, -0.1, 0.1), 
                         gamma0 = gamma0[max.spec] + runif(1, -0.1, 0.1), 
                         gamma1 = gamma1[max.spec] + runif(1, -0.1, 0.1)
)}

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1"))

MCMCconf$addSampler(target = c("gamma0", "gamma1"),
                    type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  output <- rbind(output,
                  data.frame(model = "sdm.ss",
                             parameter.type = "max.spec",
                             parameter = params[k],
                             truth = get(params[k])[target[2]],
                             mean = summary(out)[[1]][params[k],"Mean"],
                             sd = summary(out)[[1]][params[k], "SD"],
                             rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
}

#---------------------------------------#
#-Auxiliary Single Species Model Code-#
#---------------------------------------#

code <- nimbleCode({
  
  #PRIORS
  
  #Population growth
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)  
  
  #Population density
  beta0 ~ dnorm(0, 0.01) 
  beta1 ~ dnorm(0, 0.01)
  
  for(j in 1:nsites){
    
    y[1,j] ~ dpois(lambda[1,j])
    
    lambda[1,j] <- exp(beta0 + beta1 * cov.st[1,j] + log.correction)
    
  }#end j
  
  for(t in 2:nyears){
    
    log(gamma[t-1]) <- gamma0 + gamma1 * cov.t[t-1]
    
    for(j in 1:nsites){
      
      y[t,j] ~ dpois(lambda[t,j])
      
      lambda[t,j] <- exp(beta0 + log(prod(gamma[1:(t-1)])) + beta1 * cov.st[t,j] + log.correction)
      
    }#end j
    
  }#end t
  
})

#-Compile data-#

data <- list(y = y[min.spec,,],
             cov.t = cov.t,
             cov.st = cov.st,
             correction = summary(aux.out$samples2)[[1]][paste0("log.correction[", min.spec, "]"),"Mean"])

constants <- list(nyears = nyears, nsites = nsites)

#-Initial values-#

inits <- function(){list(beta0 = beta0[min.spec] + runif(1, -0.1, 0.1), 
                         beta1 = beta1[min.spec] + runif(1, -0.1, 0.1), 
                         gamma0 = gamma0[min.spec] + runif(1, -0.1, 0.1), 
                         gamma1 = gamma1[min.spec] + runif(1, -0.1, 0.1)
)}

#-Parameters to save-#

params <- c(
  "beta0",
  "beta1",
  "gamma0",
  "gamma1"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1"))

MCMCconf$addSampler(target = c("gamma0", "gamma1"),
                    type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  output <- rbind(output,
                  data.frame(model = "aux.ss",
                             parameter.type = "min.spec",
                             parameter = params[k],
                             truth = get(params[k])[target[1]],
                             mean = summary(out)[[1]][params[k],"Mean"],
                             sd = summary(out)[[1]][params[k], "SD"],
                             rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
}

#-Compile data-#

data <- list(y = y[max.spec,,],
             cov.t = cov.t,
             cov.st = cov.st,
             correction = summary(aux.out$samples2)[[1]][paste0("log.correction[", max.spec, "]"),"Mean"])

#-Initial values-#

inits <- function(){list(beta0 = beta0[max.spec] + runif(1, -0.1, 0.1), 
                         beta1 = beta1[max.spec] + runif(1, -0.1, 0.1), 
                         gamma0 = gamma0[max.spec] + runif(1, -0.1, 0.1), 
                         gamma1 = gamma1[max.spec] + runif(1, -0.1, 0.1)
)}

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1"))

MCMCconf$addSampler(target = c("gamma0", "gamma1"),
                    type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  output <- rbind(output,
                  data.frame(model = "aux.ss",
                             parameter.type = "max.spec",
                             parameter = params[k],
                             truth = get(params[k])[target[2]],
                             mean = summary(out)[[1]][params[k],"Mean"],
                             sd = summary(out)[[1]][params[k], "SD"],
                             rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
}

#---------------------------------#
#-Prior Single Species Model Code-#
#---------------------------------#

code <- nimbleCode({
  
  #PRIORS
  
  #Population growth
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)  
  
  #Population density
  beta0 ~ dnorm(0, 0.01) 
  beta1 ~ dnorm(0, 0.01)
  
  #Correction factor
  log.correction ~ dnorm(mu.log.correction, sd = sd.log.correction)
  
  for(j in 1:nsites){
    
    y[1,j] ~ dpois(lambda[1,j])
    
    lambda[1,j] <- exp(beta0 + beta1 * cov.st[1,j] + log.correction)
    
  }#end j
  
  for(t in 2:nyears){
    
    log(gamma[t-1]) <- gamma0 + gamma1 * cov.t[t-1]
    
    for(j in 1:nsites){
      
      y[t,j] ~ dpois(lambda[t,j])
      
      lambda[t,j] <- exp(beta0 + log(prod(gamma[1:(t-1)])) + beta1 * cov.st[t,j] + log.correction)
      
    }#end j
    
  }#end t
  
})

#-Compile data-#

data <- list(y = y[min.spec,,],
             cov.t = cov.t,
             cov.st = cov.st,
             mu.log.correction = summary(aux.out$samples2)[[1]][paste0("log.correction[", min.spec, "]"),"Mean"],
             sd.log.correction = summary(aux.out$samples2)[[1]][paste0("log.correction[", min.spec, "]"),"SD"])

constants <- list(nyears = nyears, nsites = nsites)

#-Initial values-#

inits <- function(){list(beta0 = beta0[min.spec] + runif(1, -0.1, 0.1), 
                         beta1 = beta1[min.spec] + runif(1, -0.1, 0.1), 
                         gamma0 = gamma0[min.spec] + runif(1, -0.1, 0.1), 
                         gamma1 = gamma1[min.spec] + runif(1, -0.1, 0.1),
                         correction = summary(aux.out$samples2)[[1]][paste0("log.correction[", min.spec, "]"),"Mean"]
)}

#-Parameters to save-#

params <- c(
  "beta0",
  "beta1",
  "gamma0",
  "gamma1",
  "log.correction"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1", "beta0", "beta1"))

MCMCconf$addSampler(target = c("gamma0", "gamma1"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("beta0", "beta1"),
                    type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  output <- rbind(output,
                  data.frame(model = "prior.ss",
                             parameter.type = "min.spec",
                             parameter = params[k],
                             truth = get(params[k])[target[1]],
                             mean = summary(out)[[1]][params[k],"Mean"],
                             sd = summary(out)[[1]][params[k], "SD"],
                             rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
}

#-Compile data-#

data <- list(y = y[max.spec,,],
             cov.t = cov.t,
             cov.st = cov.st,
             mu.correction = summary(aux.out$samples2)[[1]][paste0("log.correction[", max.spec, "]"),"Mean"],
             sd.correction = summary(aux.out$samples2)[[1]][paste0("log.correction[", max.spec, "]"),"SD"])

#-Initial values-#

inits <- function(){list(beta0 = beta0[max.spec] + runif(1, -0.1, 0.1), 
                         beta1 = beta1[max.spec] + runif(1, -0.1, 0.1), 
                         gamma0 = gamma0[max.spec] + runif(1, -0.1, 0.1), 
                         gamma1 = gamma1[max.spec] + runif(1, -0.1, 0.1),
                         correction = summary(aux.out$samples2)[[1]][paste0("log.correction[", max.spec, "]"),"Mean"]
)}

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1", "beta0", "beta1"))

MCMCconf$addSampler(target = c("gamma0", "gamma1"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("beta0", "beta1"),
                    type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  output <- rbind(output,
                  data.frame(model = "prior.ss",
                             parameter.type = "max.spec",
                             parameter = params[k],
                             truth = get(params[k])[target[2]],
                             mean = summary(out)[[1]][params[k],"Mean"],
                             sd = summary(out)[[1]][params[k], "SD"],
                             rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
}

#--------------------------------------#
#-Integrated Single Species Model Code-#
#--------------------------------------#

code <- nimbleCode({
  
  #PRIORS
  
  for(i in 1:nspecies){
    
    gamma0[i] ~ dnorm(0, 0.01)
    gamma1[i] ~ dnorm(0, 0.01)
    
    #Population density
    beta0[i] ~ dnorm(0, 0.01) 
    beta1[i] ~ dnorm(0, 0.01)
    
    #LIKELIHOOD
    
    for(j in 1:nsites){
      
      y[i,1,j] ~ dpois(lambda[i,1,j])
      
      lambda[i,1,j] <- exp(beta0[i] + beta1[i] * cov.st[1,j] + log.correction[i])
      
    }#end j
    
    for(t in 2:nyears){
      
      log(gamma[i,t-1]) <- gamma0[i] + gamma1[i] * cov.t[t-1]
      
      for(j in 1:nsites){
        
        y[i,t,j] ~ dpois(lambda[i,t,j])
        
        lambda[i,t,j] <- exp(beta0[i] + log(prod(gamma[i,1:(t-1)])) + beta1[i] * cov.st[t,j] + log.correction[i]) 
        
      }#end j
      
    }#end t
    
    log.correction[i] <- log(correction[i])
    
    correction[i] <- epsilon.hat * phi[i]/pi[i]
    
    phi.ones[i] <- 1
    
    exp.beta[i] <- exp(beta0[i])
    
    pi[i] <- exp(beta0[i])/sum(exp.beta[1:nspecies])
    
    for(k in 1:nreps){
      
      log(theta[i,k]) <- beta0[i] + beta1[i] * cov.cam[k]
      pi.star[i,k] <- theta[i,k]/THETA[k]
      phi.star[i,k] <- (phi[i] * sum(exp.beta[1:nspecies]) * exp(beta1[i] * cov.cam[k]))/THETA[k] 
      
    }
    
  }#end i
  
  for(k in 1:nreps){
    
    FF[1:nspecies,k] ~ dmulti(pi.star[1:nspecies,k], FF.total[k])
    
    #Front facing camera total abundance
    FF.total[k] ~ dpois(THETA[k])
    
    #Community-wide expected abundance
    THETA[k] <- sum(theta[1:nspecies,k])
    
    for(o in 1:nobs){
      
      OBS[1:nspecies,k,o] ~ dmulti(phi.star[1:nspecies,k], OBS.total[k,o])
      
      OBS.total[k,o] ~ dpois(THETA[k] * epsilon.hat)
      
    }#end o
    
  }#end k
  
  #Composition of latent abundance (corrected for imperfect detection)
  phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  
  #Derived product of movement and detection
  epsilon.hat ~ dnorm(0, 0.01)
  
})

#-Compile data-#

data <- list(y = y,
             cov.t = cov.t,
             cov.st = cov.st,
             cov.cam = cov.cam,
             FF = FF,
             FF.total = apply(FF, 2, sum),
             OBS = OBS,
             OBS.total = apply(OBS, c(2,3), sum))

constants <- list(nspecies = nspecies, nyears = nyears, nsites = nsites,
                  nreps = nreps, nobs = 2)

#-Initial values-#

inits <- function(){list(beta0 = beta0 + runif(nspecies, -0.1, 0.1), 
                         beta1 = beta1 + runif(nspecies, -0.1, 0.1), 
                         gamma0 = gamma0 + runif(nspecies, -0.1, 0.1), 
                         gamma1 = gamma1 + runif(nspecies, -0.1, 0.1),
                         epsilon.hat = alpha*p  + runif(1, -0.1, 0.1)
)}

#-Parameters to save-#

params <- c(
  "beta0",
  "beta1",
  "gamma0",
  "gamma1",
  "epsilon.hat",
  "pi",
  "phi",
  "log.correction"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1", "beta0", "epsilon.hat"))

MCMCconf$addSampler(target = c("beta0", "epsilon.hat"),
                    type = "AF_slice")

for(i in 1:nspecies){
  MCMCconf$addSampler(target = c(paste0("gamma0[", i, "]"), paste0("gamma1[", i, "]")),
                      type = "AF_slice")
}

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  if(length(get(params[k])) == 1){
    output <- rbind(output,
                    data.frame(model = "integrated.ss",
                               parameter.type = "community",
                               parameter = params[k],
                               truth = get(params[k]),
                               mean = summary(out)[[1]][params[k],"Mean"],
                               sd = summary(out)[[1]][params[k], "SD"],
                               rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
  }else{
    for(j in 1:2){
      output <- rbind(output,
                      data.frame(model = "integrated.ss",
                                 parameter.type = c("min.spec", "max.spec")[j],
                                 parameter = params[k],
                                 truth = get(params[k])[target[j]],
                                 mean = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"Mean"],
                                 sd = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"SD"],
                                 rhat = as.numeric(coda::gelman.diag(out[1:nc][,paste0(params[k], "[", target[j], "]")])[[1]][,1])))
    }
  }
}

#----------------------------------#
#-Single Data Community Model Code-#
#----------------------------------#

code <- nimbleCode({
  
  #PRIORS
  
  #Population growth
  mu.g0 ~ dnorm(0, 0.01)
  tau.g0 ~ dgamma(1, 0.001)
  sd.g0 <- 1/sqrt(tau.g0)
  
  mu.g1 ~ dnorm(0, 0.01)
  tau.g1 ~ dgamma(1, 0.001)
  sd.g1 <- 1/sqrt(tau.g1)
  
  #Population density
  mu.b0 ~ dnorm(0, 0.01)
  tau.b0 ~ dgamma(1, 0.001)
  sd.b0 <- 1/sqrt(tau.b0)
  
  mu.b1 ~ dnorm(0, 0.01)
  tau.b1 ~ dgamma(1, 0.001)
  sd.b1 <- 1/sqrt(tau.b1)
  
  for(i in 1:nspecies){
    
    gamma0[i] ~ dnorm(mu.g0, tau.g0)
    gamma1[i] ~ dnorm(mu.g1, tau.g1)
    
    #Population density
    beta0[i] ~ dnorm(mu.b0, tau.b0) 
    beta1[i] ~ dnorm(mu.b1, tau.b1)
    
    #LIKELIHOOD
    
    for(j in 1:nsites){
      
      y[i,1,j] ~ dpois(lambda[i,1,j])
      
      lambda[i,1,j] <- exp(beta0[i] + beta1[i] * cov.st[1,j])
      
    }#end j
    
    for(t in 2:nyears){
      
      log(gamma[i,t-1]) <- gamma0[i] + gamma1[i] * cov.t[t-1]
      
      for(j in 1:nsites){
        
        y[i,t,j] ~ dpois(lambda[i,t,j])
        
        lambda[i,t,j] <- exp(beta0[i] + log(prod(gamma[i,1:(t-1)])) + beta1[i] * cov.st[t,j])
        
      }#end j
      
    }#end t
    
  }#end i
  
})

#-Compile data-#

data <- list(y = y,
             cov.t = cov.t,
             cov.st = cov.st)

constants <- list(nspecies = nspecies, nyears = nyears, nsites = nsites)

#-Initial values-#

inits <- function(){list(beta0 = beta0 + runif(nspecies, -0.1, 0.1), 
                         mu.b0 = mu.b0 + runif(1, -0.1, 0.1), 
                         sd.b0 = sd.b0 + runif(1, 0, 0.25),
                         beta1 = beta1 + runif(nspecies, -0.1, 0.1), 
                         mu.b1 = mu.b1 + runif(1, -0.1, 0.1), 
                         sd.b1 = sd.b1 + runif(1, -0.1, 0.1),
                         gamma0 = gamma0 + runif(nspecies, -0.1, 0.1), 
                         mu.g0 = mu.g0 + runif(1, -0.1, 0.1), 
                         sd.g0 = sd.g0 + runif(1, -0.1, 0.1),
                         gamma1 = gamma1 + runif(nspecies, -0.1, 0.1),
                         mu.g1 = mu.g0 + runif(1, -0.1, 0.1), 
                         sd.g1 = sd.g1 + runif(1, -0.1, 0.1)
)}

#-Parameters to save-#

params <- c(
  "mu.b0",
  "sd.b0",
  "mu.b1",
  "sd.b1",
  "beta0",
  "beta1",
  "mu.g0",
  "sd.g0",
  "mu.g1",
  "sd.g1",
  "gamma0",
  "gamma1"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1"))

for(i in 1:nspecies){
  MCMCconf$addSampler(target = c(paste0("gamma0[", i, "]"), paste0("gamma1[", i, "]")),
                      type = "AF_slice")
}
MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

for(k in 1:length(params)){
  if(length(get(params[k])) == 1){
    output <- rbind(output,
                    data.frame(model = "sdm.cm",
                               parameter.type = "community",
                               parameter = params[k],
                               truth = get(params[k]),
                               mean = summary(out)[[1]][params[k],"Mean"],
                               sd = summary(out)[[1]][params[k], "SD"],
                               rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
  }else{
    for(j in 1:2){
      output <- rbind(output,
                      data.frame(model = "sdm.cm",
                                 parameter.type = c("min.spec", "max.spec")[j],
                                 parameter = params[k],
                                 truth = get(params[k])[target[j]],
                                 mean = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"Mean"],
                                 sd = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"SD"],
                                 rhat = as.numeric(coda::gelman.diag(out[1:nc][,paste0(params[k], "[", target[j], "]")])[[1]][,1])))
    }
  }
}

#-------------------------------#
#-Aux Data Community Model Code-#
#-------------------------------#

code <- nimbleCode({
  
  #PRIORS
  
  #Population growth
  mu.g0 ~ dnorm(0, 0.01)
  tau.g0 ~ dgamma(1, 0.001)
  sd.g0 <- 1/sqrt(tau.g0)
  
  mu.g1 ~ dnorm(0, 0.01)
  tau.g1 ~ dgamma(1, 0.001)
  sd.g1 <- 1/sqrt(tau.g1)
  
  #Population density
  mu.b0 ~ dnorm(0, 0.01)
  tau.b0 ~ dgamma(1, 0.001)
  sd.b0 <- 1/sqrt(tau.b0)
  
  mu.b1 ~ dnorm(0, 0.01)
  tau.b1 ~ dgamma(1, 0.001)
  sd.b1 <- 1/sqrt(tau.b1)
  
  for(i in 1:nspecies){
    
    gamma0[i] ~ dnorm(mu.g0, tau.g0)
    gamma1[i] ~ dnorm(mu.g1, tau.g1)
    
    #Population density
    beta0[i] ~ dnorm(mu.b0, tau.b0) 
    beta1[i] ~ dnorm(mu.b1, tau.b1)
    
    #LIKELIHOOD
    
    for(j in 1:nsites){
      
      y[i,1,j] ~ dpois(lambda[i,1,j])
      
      lambda[i,1,j] <- exp(beta0[i] + beta1[i] * cov.st[1,j] + log.correction[i])
      
    }#end j
    
    for(t in 2:nyears){
      
      log(gamma[i,t-1]) <- gamma0[i] + gamma1[i] * cov.t[t-1]
      
      for(j in 1:nsites){
        
        y[i,t,j] ~ dpois(lambda[i,t,j])
        
        lambda[i,t,j] <- exp(beta0[i] + log(prod(gamma[i,1:(t-1)])) + beta1[i] * cov.st[t,j] + log.correction[i])
        
      }#end j
      
    }#end t
    
  }#end i
  
})

#-Compile data-#

data <- list(y = y,
             cov.t = cov.t,
             cov.st = cov.st,
             log.correction = summary(aux.out$samples2)[[1]][,"Mean"])

constants <- list(nspecies = nspecies, nyears = nyears, nsites = nsites)

#-Initial values-#

inits <- function(){list(beta0 = beta0 + runif(nspecies, -0.1, 0.1), 
                         mu.b0 = mu.b0 + runif(1, -0.1, 0.1), 
                         sd.b0 = sd.b0 + runif(1, 0, 0.25),
                         beta1 = beta1 + runif(nspecies, -0.1, 0.1), 
                         mu.b1 = mu.b1 + runif(1, -0.1, 0.1), 
                         sd.b1 = sd.b1 + runif(1, -0.1, 0.1),
                         gamma0 = gamma0 + runif(nspecies, -0.1, 0.1), 
                         mu.g0 = mu.g0 + runif(1, -0.1, 0.1), 
                         sd.g0 = sd.g0 + runif(1, -0.1, 0.1),
                         gamma1 = gamma1 + runif(nspecies, -0.1, 0.1),
                         mu.g1 = mu.g0 + runif(1, -0.1, 0.1), 
                         sd.g1 = sd.g1 + runif(1, -0.1, 0.1)
)}

#-Parameters to save-#

params <- c(
  "mu.b0",
  "sd.b0",
  "mu.b1",
  "sd.b1",
  "beta0",
  "beta1",
  "mu.g0",
  "sd.g0",
  "mu.g1",
  "sd.g1",
  "gamma0",
  "gamma1"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1"))

for(i in 1:nspecies){
  MCMCconf$addSampler(target = c(paste0("gamma0[", i, "]"), paste0("gamma1[", i, "]")),
                      type = "AF_slice")
}

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  if(length(get(params[k])) == 1){
    output <- rbind(output,
                    data.frame(model = "aux.cm",
                               parameter.type = "community",
                               parameter = params[k],
                               truth = get(params[k]),
                               mean = summary(out)[[1]][params[k],"Mean"],
                               sd = summary(out)[[1]][params[k], "SD"],
                               rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
  }else{
    for(j in 1:2){
      output <- rbind(output,
                      data.frame(model = "aux.cm",
                                 parameter.type = c("min.spec", "max.spec")[j],
                                 parameter = params[k],
                                 truth = get(params[k])[target[j]],
                                 mean = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"Mean"],
                                 sd = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"SD"],
                                 rhat = as.numeric(coda::gelman.diag(out[1:nc][,paste0(params[k], "[", target[j], "]")])[[1]][,1])))
    }
  }
}

#---------------------------------#
#-Prior Data Community Model Code-#
#---------------------------------#

code <- nimbleCode({
  
  #PRIORS
  
  #Population growth
  mu.g0 ~ dnorm(0, 0.01)
  tau.g0 ~ dgamma(1, 0.001)
  sd.g0 <- 1/sqrt(tau.g0)
  
  mu.g1 ~ dnorm(0, 0.01)
  tau.g1 ~ dgamma(1, 0.001)
  sd.g1 <- 1/sqrt(tau.g1)
  
  #Population density
  mu.b0 ~ dnorm(0, 0.01)
  tau.b0 ~ dgamma(1, 0.001)
  sd.b0 <- 1/sqrt(tau.b0)
  
  mu.b1 ~ dnorm(0, 0.01)
  tau.b1 ~ dgamma(1, 0.001)
  sd.b1 <- 1/sqrt(tau.b1)
  
  for(i in 1:nspecies){
    
    #Correction factor
    log.correction[i] ~ dnorm(mu.log.correction[i], sd = sd.log.correction[i])
    
    #Population growth
    gamma0[i] ~ dnorm(mu.g0, tau.g0)
    gamma1[i] ~ dnorm(mu.g1, tau.g1)
    
    #Population density
    beta0[i] ~ dnorm(mu.b0, tau.b0) 
    beta1[i] ~ dnorm(mu.b1, tau.b1)
    
    #LIKELIHOOD
    
    for(j in 1:nsites){
      
      y[i,1,j] ~ dpois(lambda[i,1,j])
      
      lambda[i,1,j] <- exp(beta0[i] + beta1[i] * cov.st[1,j] + log.correction[i])
      
    }#end j
    
    for(t in 2:nyears){
      
      log(gamma[i,t-1]) <- gamma0[i] + gamma1[i] * cov.t[t-1]
      
      for(j in 1:nsites){
        
        y[i,t,j] ~ dpois(lambda[i,t,j])
        
        lambda[i,t,j] <- exp(beta0[i] + log(prod(gamma[i,1:(t-1)])) + beta1[i] * cov.st[t,j] + log.correction[i])
        
      }#end j
      
    }#end t
    
  }#end i
  
})

#-Compile data-#

data <- list(y = y,
             cov.t = cov.t,
             cov.st = cov.st,
             mu.log.correction = summary(aux.out$samples2)[[1]][,"Mean"],
             sd.log.correction = summary(aux.out$samples2)[[1]][,"SD"])

constants <- list(nspecies = nspecies, nyears = nyears, nsites = nsites)

#-Initial values-#

inits <- function(){list(beta0 = beta0 + runif(nspecies, -0.1, 0.1), 
                         mu.b0 = mu.b0 + runif(1, -0.1, 0.1), 
                         sd.b0 = sd.b0 + runif(1, 0, 0.25),
                         beta1 = beta1 + runif(nspecies, -0.1, 0.1), 
                         mu.b1 = mu.b1 + runif(1, -0.1, 0.1), 
                         sd.b1 = sd.b1 + runif(1, -0.1, 0.1),
                         gamma0 = gamma0 + runif(nspecies, -0.1, 0.1), 
                         mu.g0 = mu.g0 + runif(1, -0.1, 0.1), 
                         sd.g0 = sd.g0 + runif(1, -0.1, 0.1),
                         gamma1 = gamma1 + runif(nspecies, -0.1, 0.1),
                         mu.g1 = mu.g0 + runif(1, -0.1, 0.1), 
                         sd.g1 = sd.g1 + runif(1, -0.1, 0.1),
                         log.correction = summary(aux.out$samples2)[[1]][,"Mean"]
)}

#-Parameters to save-#

params <- c(
  "mu.b0",
  "sd.b0",
  "mu.b1",
  "sd.b1",
  "beta0",
  "beta1",
  "mu.g0",
  "sd.g0",
  "mu.g1",
  "sd.g1",
  "gamma0",
  "gamma1",
  "log.correction"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1", "beta0", "log.correction"))

for(i in 1:nspecies){
  MCMCconf$addSampler(target = c(paste0("gamma0[", i, "]"), paste0("gamma1[", i, "]")),
                      type = "AF_slice")
  MCMCconf$addSampler(target = c(paste0("beta0[", i, "]"), paste0("log.correction[", i, "]")),
                      type = "AF_slice")
}

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  if(length(get(params[k])) == 1){
    output <- rbind(output,
                    data.frame(model = "prior.cm",
                               parameter.type = "community",
                               parameter = params[k],
                               truth = get(params[k]),
                               mean = summary(out)[[1]][params[k],"Mean"],
                               sd = summary(out)[[1]][params[k], "SD"],
                               rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
  }else{
    for(j in 1:2){
      output <- rbind(output,
                      data.frame(model = "prior.cm",
                                 parameter.type = c("min.spec", "max.spec")[j],
                                 parameter = params[k],
                                 truth = get(params[k])[target[j]],
                                 mean = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"Mean"],
                                 sd = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"SD"],
                                 rhat = as.numeric(coda::gelman.diag(out[1:nc][,paste0(params[k], "[", target[j], "]")])[[1]][,1])))
    }
  }
}

# if(any(output %>% filter(model == "prior.cm") %>% select(rhat) %>%.$rhat > 1.1)){
#   nonconverged <- list(
#     out = out,
#     output = output %>% filter(model == "prior.cm")
#   )
#   nc.ID <- length(list.files("./nonconverged/priorcm/")) + 1
#   save(nonconverged, file = paste("./nonconverged/priorcm/nc", nc.ID, ".Rds", sep=""))
# }

#---------------------------------#
#-Integrated Community Model Code-#
#---------------------------------#

code <- nimbleCode({
  
  #PRIORS
  
  #Population growth
  mu.g0 ~ dnorm(0, 0.01)
  tau.g0 ~ dgamma(1, 0.001)
  sd.g0 <- 1/sqrt(tau.g0)
  
  mu.g1 ~ dnorm(0, 0.01)
  tau.g1 ~ dgamma(1, 0.001)
  sd.g1 <- 1/sqrt(tau.g1)
  
  #Population density
  mu.b0 ~ dnorm(0, 0.01)
  tau.b0 ~ dgamma(1, 0.001)
  sd.b0 <- 1/sqrt(tau.b0)
  
  mu.b1 ~ dnorm(0, 0.01)
  tau.b1 ~ dgamma(1, 0.001)
  sd.b1 <- 1/sqrt(tau.b1)
  
  for(i in 1:nspecies){
    
    gamma0[i] ~ dnorm(mu.g0, tau.g0)
    gamma1[i] ~ dnorm(mu.g1, tau.g1)
    
    #Population density
    beta0[i] ~ dnorm(mu.b0, tau.b0) 
    beta1[i] ~ dnorm(mu.b1, tau.b1)
    
    #LIKELIHOOD
    
    for(j in 1:nsites){
      
      y[i,1,j] ~ dpois(lambda[i,1,j])
      
      lambda[i,1,j] <- exp(beta0[i] + beta1[i] * cov.st[1,j] + log.correction[i])
      
    }#end j
    
    for(t in 2:nyears){
      
      log(gamma[i,t-1]) <- gamma0[i] + gamma1[i] * cov.t[t-1]
      
      for(j in 1:nsites){
        
        y[i,t,j] ~ dpois(lambda[i,t,j])
        
        lambda[i,t,j] <- exp(beta0[i] + log(prod(gamma[i,1:(t-1)])) + beta1[i] * cov.st[t,j] + log.correction[i])
        
      }#end j
      
    }#end t
    
    log.correction[i] <- log(correction[i])
    
    correction[i] <- epsilon.hat * phi[i]/pi[i]
    
    phi.ones[i] <- 1
    
    exp.beta[i] <- exp(beta0[i])
    
    pi[i] <- exp(beta0[i])/sum(exp.beta[1:nspecies])
    
    for(k in 1:nreps){
      
      log(theta[i,k]) <- beta0[i] + beta1[i] * cov.cam[k]
      pi.star[i,k] <- theta[i,k]/THETA[k]
      phi.star[i,k] <- (phi[i] * sum(exp.beta[1:nspecies]) * exp(beta1[i] * cov.cam[k]))/THETA[k] 
      
    }
    
  }#end i
  
  for(k in 1:nreps){
    
    FF[1:nspecies,k] ~ dmulti(pi.star[1:nspecies,k], FF.total[k])
    
    #Front facing camera total abundance
    FF.total[k] ~ dpois(THETA[k])
    
    #Community-wide expected abundance
    THETA[k] <- sum(theta[1:nspecies,k])
    
    for(o in 1:nobs){
      
      OBS[1:nspecies,k,o] ~ dmulti(phi.star[1:nspecies,k], OBS.total[k,o])
      
      OBS.total[k,o] ~ dpois(THETA[k] * epsilon.hat)
      
    }#end o
    
  }#end k
  
  #Composition of latent abundance (corrected for imperfect detection)
  phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  
  #Derived product of movement and detection
  epsilon.hat ~ dnorm(0, 0.01)
  
})

#-Compile data-#

data <- list(y = y,
             cov.t = cov.t,
             cov.st = cov.st,
             cov.cam = cov.cam,
             FF = FF,
             FF.total = apply(FF, 2, sum),
             OBS = OBS,
             OBS.total = apply(OBS, c(2,3), sum))

constants <- list(nspecies = nspecies, nyears = nyears, nsites = nsites,
                  nreps = nreps, nobs = 2)

#-Initial values-#

inits <- function(){list(beta0 = beta0 + runif(nspecies, -0.1, 0.1), 
                         mu.b0 = mu.b0 + runif(1, -0.1, 0.1), 
                         sd.b0 = sd.b0 + runif(1, 0, 0.25),
                         beta1 = beta1 + runif(nspecies, -0.1, 0.1), 
                         mu.b1 = mu.b1 + runif(1, -0.1, 0.1), 
                         sd.b1 = sd.b1 + runif(1, -0.1, 0.1),
                         gamma0 = gamma0 + runif(nspecies, -0.1, 0.1), 
                         mu.g0 = mu.g0 + runif(1, -0.1, 0.1), 
                         sd.g0 = sd.g0 + runif(1, -0.1, 0.1),
                         gamma1 = gamma1 + runif(nspecies, -0.1, 0.1),
                         mu.g1 = mu.g0 + runif(1, -0.1, 0.1), 
                         sd.g1 = sd.g1 + runif(1, -0.1, 0.1),
                         epsilon.hat = alpha*p  + runif(1, -0.1, 0.1)
)}

#-Parameters to save-#

params <- c(
  "mu.b0",
  "sd.b0",
  "mu.b1",
  "sd.b1",
  "beta0",
  "beta1",
  "mu.g0",
  "sd.g0",
  "mu.g1",
  "sd.g1",
  "gamma0",
  "gamma1",
  "epsilon.hat",
  "pi",
  "phi",
  "log.correction"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMCconf$removeSampler(target = c("gamma0", "gamma1", "beta0", "epsilon.hat"))

MCMCconf$addSampler(target = c("beta0", "epsilon.hat"),
                    type = "AF_slice")

for(i in 1:nspecies){
  MCMCconf$addSampler(target = c(paste0("gamma0[", i, "]"), paste0("gamma1[", i, "]")),
                      type = "AF_slice")
}

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 12000
nb <- 2000
nt <- 10

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Output-#

for(k in 1:length(params)){
  if(length(get(params[k])) == 1){
    output <- rbind(output,
                    data.frame(model = "integrated.cm",
                               parameter.type = "community",
                               parameter = params[k],
                               truth = get(params[k]),
                               mean = summary(out)[[1]][params[k],"Mean"],
                               sd = summary(out)[[1]][params[k], "SD"],
                               rhat = as.numeric(coda::gelman.diag(out[1:nc][,params[k]])[[1]][,1])))
  }else{
    for(j in 1:2){
      output <- rbind(output,
                      data.frame(model = "integrated.cm",
                                 parameter.type = c("min.spec", "max.spec")[j],
                                 parameter = params[k],
                                 truth = get(params[k])[target[j]],
                                 mean = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"Mean"],
                                 sd = summary(out)[[1]][paste0(params[k], "[", target[j], "]"),"SD"],
                                 rhat = as.numeric(coda::gelman.diag(out[1:nc][,paste0(params[k], "[", target[j], "]")])[[1]][,1])))
    }
  }
}


# if(any(output %>% filter(model == "integrated.cm") %>% select(rhat) %>%.$rhat > 1.1)){
#   nonconverged <- list(
#     out = out,
#     output = output %>% filter(model == "integrated.cm")
#   )
#   nc.ID <- length(list.files("./nonconverged/icm/")) + 1
#   save(nonconverged, file = paste("./nonconverged/icm/nc", nc.ID, ".Rds", sep=""))
# }

#-Export output-#
ID <- length(list.files("./PostAnalysis/simulation_output/")) + 1
save(output, file = paste(".PostAnalysis/simulation_output/output", ID, ".Rds", sep=""))

#}
