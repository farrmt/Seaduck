library(nimble)
library(coda)

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
  phi.pi[sample(x = 1:nspecies^2, size = nmissID)] <- runif(nmissID, 0, 0.25)
  diag(phi.pi) <- runif(nspecies, 0.8, 1)
  phi.pi <- phi.pi/apply(phi.pi, MARGIN = 1, sum)
  return(phi.pi)
}

#Generate composition proportions
comp.fun <- function(nspecies)
{
  pi <- runif(nspecies, 0, 1)
  pi <- pi/sum(pi)
  return(pi)
}

#Number of years
nyears <- 10

#Number of sites
nsites <- 500

#Number of species
nspecies <- 5

#Expected population growth rate intercept
gamma0 <- log(runif(nspecies, 0.9, 1.1))

#Expected population growth rate covariate effect
gamma1 <- runif(nspecies, -0.25, 0.25)

#Expected population growth rate covariate
cov.t <- rnorm(nyears - 1, mean = 0, sd = 1)

#Expected population growth rate
gamma <- matrix(data = NA, nrow = nspecies, ncol = nyears - 1)

for(i in 1:nspecies){
  gamma[i,] <- exp(gamma0[i] + gamma1[i] * cov.t)
}

#Expected population abundance
LAMBDA <- matrix(data = NA, nrow = nspecies, ncol = nyears)

#Expected population abundance in year 1
LAMBDA[,1] <- runif(nspecies, 10000, 20000)

#Expected population abundance in remaining years
for(t in 2:nyears){
  LAMBDA[,t] <- LAMBDA[,t-1] * gamma[,t-1]
}

#Latent population abundance
N <- matrix(data = NA, nrow = nspecies, ncol = nyears)

for(i in 1:nspecies){
  N[i,] <- rpois(nyears, LAMBDA[i,])
}

#Occupancy/use
omega0 <- runif(1, 0.1, 0.5)

#Occupancy/use covariate
cov.use <- rnorm(nsites, mean = 0, sd = 1)

#Occupancy/use covariate effect
omega1 <- runif(1, -1, 1)

#Expected population abudance at the site-level
omega <- expit(logit(omega0) + omega1 * cov.use)

#Latent use
z <- rbinom(n = nsites, size = 1, prob = omega)

#Expected population abundance at the site-level
lambda <- array(NA, dim = c(nspecies, nyears, nsites))

#Expected population abundance intercept at the site-level
beta0 <- log(LAMBDA[,1]/sum(z))
beta0.org <- beta0

#Spatiotemporal covariate effects
beta1 <- runif(nspecies, -1, 1)

#Spatiotemporal covariate
cov.st <- matrix(NA, nrow = nyears, ncol = nsites)

for(t in 1:nyears){
  cov.st[t,1:nsites] <- rnorm(nsites, mean = 0, sd = 1)
  cov.st[t,1:nsites] <- scale(cov.st[t,1:nsites])
}

#Proportion of total population abundance at a site
site.prob <- array(NA, dim = c(nspecies, nyears, nsites))

#Latent site-level population abundance
n <- array(NA, dim = c(nspecies, nyears, nsites))

#Generate values for the above
for(i in 1:nspecies){
  lambda[i,1,1:nsites] <- z[1:nsites] * exp(beta0[i] + beta1[i] * cov.st[1,1:nsites])
  lambda[i,1,1:nsites] <- lambda[i,1,1:nsites]/sum(lambda[i,1,1:nsites]) * LAMBDA[i,1]
  site.prob[i,1,1:nsites] <- lambda[i,1,1:nsites]/LAMBDA[i,1]
  n[i,1,1:nsites] <- rmultinom(1, N[i,1], site.prob[i,1,1:nsites])
  for(t in 1:nyears){
    lambda[i,t,1:nsites] <- z[1:nsites] * exp(beta0[i] + beta1[i] * cov.st[t,1:nsites] + log(prod(gamma[i,1:(t-1)])))
    lambda[i,t,1:nsites] <- lambda[i,t,1:nsites]/sum(lambda[i,t,1:nsites]) * LAMBDA[i,t]
    site.prob[i,t,1:nsites] <- lambda[i,t,1:nsites]/LAMBDA[i,t]
    n[i,t,1:nsites] <- rmultinom(1, N[i,t], site.prob[i,t,1:nsites])
  }
}

for(i in 1:nspecies){
  beta0[i] <- mean(log(lambda[i,1,z==1]) - beta1[i] * cov.st[1,z==1])
}

#-Observation process-#

#Detection probability
p <- runif(1, 0, 1)

#Movement rate
alpha <- runif(1, 0.5, 1)

#Expected miss ID rate
phi.pi <- missID.fun(nspecies, nmissID = 4)

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
nreps <- 100

#Community composition
pi <- comp.fun(nspecies)

#Community epected abundance
THETA <- runif(1, 300, 1000)

#Species expected abundance
theta <- THETA * pi

#Front facing and point of view camera data
FF <- array(NA, dim = c(nspecies, nreps))

#Observer data
OBS <- array(NA, dim = c(nspecies, nreps, 2))

for(k in 1:nreps){
  FF[,k] <- rpois(nspecies, theta) 
  tmp <- NULL 
  for(i in 1:nspecies){
    tmp <- cbind(tmp, rmultinom(1, FF[i,k], phi.pi[i,]))
  }#end i
  OBS[,k,1] <- rbinom(n = nspecies, apply(tmp, MARGIN = 1, sum), prob = p * alpha)
  OBS[,k,2] <- rbinom(n = nspecies, apply(tmp, MARGIN = 1, sum), prob = p * alpha)
}#end k


#-Nimble Code-#

code <- nimbleCode({
  
  #PRIORS
  
  omega0 ~ dbeta(1, 1)
  
  omega1~ dnorm(0, 0.01)
  
  for(j in 1:nsites){
    z[j] ~ dbern(omega[j])
    
    logit(omega[j]) <- logit(omega0) + omega1 * cov.use[j]
  }
  
  for(i in 1:nspecies){
    
    gamma1[i] ~ dnorm(0, 0.01)
    
    gamma0[i] ~ dnorm(0, 0.01)
    
    #beta0[i,1] ~ dnorm(0, 0.01)
    beta0[i] ~ dnorm(0, 0.01)
    
    beta1[i] ~ dnorm(0, 0.01)
    
    #LIKELIHOOD

    Y[i,1] ~ dpois(LAMBDA[i,1] * correction[i])
    
    LAMBDA[i,1] <- sum(lambda[i,1,1:nsites])
    
    for(j in 1:nsites){
      
      y[i,1,j] ~ dbin(site.prob[i,1,j], Y[i,1])
      
      site.prob[i,1,j] <- lambda[i,1,j]/LAMBDA[i,1]
      
      lambda[i,1,j] <- z[j] * exp(beta0[i] + beta1[i] * cov.st[1,j])
      
    }#end j
      
    for(t in 2:nyears){
      
      log(gamma[i,t-1]) <- gamma0[i] + gamma1[i] * cov.t[t-1]
      
      Y[i,t] ~ dpois(LAMBDA[i,t] * correction[i])
      
      LAMBDA[i,t] <- LAMBDA[i,t-1] * gamma[i,t-1]
      
      for(j in 1:nsites){
        
        y[i,t,j] ~ dbin(site.prob[i,t,j], Y[i,t])
        
        site.prob[i,t,j] <- lambda[i,t,j]/LAMBDA[i,t]
        
        lambda[i,t,j] <- z[j] * exp(beta0[i] + log(prod(gamma[i,1:(t-1)])) + beta1[i] * cov.st[t,j])
        
      }#end j
        
    }#end t
    
    pi[i] <- theta[i]/THETA
    phi.ones[i] <- 1
    
    log(theta[i]) <- theta0[i]
    theta0[i] ~ dnorm(0, 0.01)
    
    correction[i] <- epsilon.hat * phi[i]/pi[i]
      
  }#end i
  
  for(k in 1:nreps){
    
    FF[1:nspecies,k] ~ dmulti(pi[1:nspecies], FF.total[k])
    
    #Front facing camera total abundance
    FF.total[k] ~ dpois(THETA)
    
    for(o in 1:nobs){
      
      OBS[1:nspecies,k,o] ~ dmulti(phi[1:nspecies], OBS.total[k,o])
      
      OBS.total[k,o] ~ dpois(THETA * epsilon.hat)
      
    }#end o
    
  }#end k
  
  #Composition of latent abundance (corrected for imperfect detection)
  phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  
  #Derived product of movement and detection
  epsilon.hat ~ dnorm(0, 0.01)
  
  #Community-wide expected abundance
  THETA <- sum(theta[1:nspecies])

})
      
#-Compile data-#

data <- list(y = y,
             Y = apply(y, c(1,2), sum),
             cov.t = cov.t,
             cov.st = cov.st,
             cov.use = cov.use,
             FF = FF,
             FF.total = apply(FF, 2, sum),
             OBS = OBS,
             OBS.total = apply(OBS, c(2,3), sum))

constants <- list(nspecies = nspecies, nyears = nyears, nsites = nsites,
                  nreps = nreps, nobs = 2)

#-Initial values-#

inits <- function(){list(
                         #n = apply(y + 1, MARGIN = c(1,2,3), max),
                         #N = apply(apply(y + 1, MARGIN = c(1,2,4), sum), MARGIN = c(1,2), max),
                         #p = runif(1, p - 0.1*p, p + 0.1*p),
                         beta0 = beta0,
                         beta1 = beta1,
                         gamma0 = gamma0,
                         gamma1 = gamma1,
                         omega0 = omega0,
                         omega1 = omega1,
                         z = z,
                         epsilon.hat = alpha*p,
                         pi = pi
                         )}

#-Parameters to save-#

params <- c(
  "beta0",
  "beta1",
  "gamma0",
  "gamma1",
  "omega0", 
  "omega1",
  "epsilon.hat",
  "pi",
  "phi"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 20000
nb <- 10000
nt <- 5

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)
