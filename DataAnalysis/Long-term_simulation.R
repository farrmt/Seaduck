library(nimble)

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

#Number of years
nyears <- 10

#Number of sites
nsites <- 1000

#Number of species
nspecies <- 4

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
beta <- lambda <- array(NA, dim = c(nspecies, nyears, nsites))

#Expected population abundance intercept at the site-level
beta0 <- log(LAMBDA/nsites)
beta0.org <- beta0

#Spatiotemporal covariate effects
beta1 <- runif(nspecies, -1, 1)
# beta1 <- rep(runif(1, -1, 1), nspecies)

#Spatiotemporal covariate
cov.st <- matrix(NA, nrow = nyears, ncol = nsites)

for(t in 1:nyears){
  cov.st[t,1:nsites] <- rnorm(nsites, mean = 0, sd = 1)
  cov.st[t,1:nsites] <- scale(cov.st[t,1:nsites])
  # cov.st[t,1:nsites] <- rlnorm(nsites, meanlog = 0, sdlog = 1)
}

#Proportion of total population abundance at a site
pi <- array(NA, dim = c(nspecies, nyears, nsites))

#Latent site-level population abundance
n <- array(NA, dim = c(nspecies, nyears, nsites))

#Generate values for the above
for(i in 1:nspecies){
  for(t in 1:nyears){
    beta[i,t,1:nsites] <- exp(beta0[i,t] + beta1[i] * cov.st[t,1:nsites])
    lambda[i,t,1:nsites] <- beta[i,t,1:nsites] * z[1:nsites]
    lambda[i,t,1:nsites] <- lambda[i,t,1:nsites]/sum(lambda[i,t,1:nsites]) * LAMBDA[i,t]
    pi[i,t,1:nsites] <- lambda[i,t,1:nsites]/LAMBDA[i,t]
    n[i,t,1:nsites] <- rmultinom(1, N[i,t], pi[i,t,1:nsites])
    beta0[i,t] <- mean(log(lambda[i,t,z==1]) - beta1[i] * cov.st[t,z==1])
  }
}


#-Observation process-#

#nreps <- 5

# #Miss ID
# psi <- 
# 
#Detection probability
#p <- runif(1, 0.5, 1)

# #Movement rate
# movement <- runif(1, 0.5, 1)
# 
#Observed counts
# y <- array(NA, dim = c(nspecies, nyears, nsites, nreps))
# 
# for(i in 1:nspecies){
#   for(t in 1:nyears){
#     for(j in 1:nsites){
#       for(k in 1:nreps){
#         y[i,t,j,k] <- rbinom(1, n[i,t,j], p)
#       }#end k
# 
# #       n.avail[i,t,j] <- n[i,t,j] * movement
# #       
# #       confusion.matrix <- NULL
# #       
# #       for(k in 1:nspecies){
# #         
# #         confusion.matrix <- cbind(confusion.matrix, rmultinom(1, n.avail[i,t,j], phi.psi[i,]))
# #         
# #       }
# #       
# #       n.missID[i,t,j] <- apply(confusion.matrix, MARGIN = 1, sum)
# #       
#       #y[i,t,j] <- rbinom(1, n[i,t,j], p)
# 
#     }#end j
#   }#end t
# }#end i


#-Nimble Code-#

code <- nimbleCode({
  
  #PRIORS
  
  #p ~ dunif(0, 1)
  #p ~ dbeta(1, 1)
  omega0 ~ dbeta(1, 1)
  
  omega1~ dnorm(0, 0.01)
  
  for(j in 1:nsites){
    z[j] ~ dbern(omega[j])
    
    logit(omega[j]) <- logit(omega0) + omega1 * cov.use[j]
  }
  
  for(i in 1:nspecies){
    
    gamma1[i] ~ dnorm(0, 0.01)
    
    gamma0[i] ~ dnorm(0, 0.01)
    
    beta0[i,1] ~ dnorm(0, 0.01)
    
    beta1[i] ~ dnorm(0, 0.01)
    
    #LIKELIHOOD
    # 
    # for(k in 1:nreps){
    #   Y[i,1,k] ~ dbin(p, N[i,1])
    # }#end k

    N[i,1] ~ dpois(LAMBDA[i,1])
    
    LAMBDA[i,1] <- sum(lambda[i,1,1:nsites])
    
    for(j in 1:nsites){
      
      n[i,1,j] ~ dbin(pi[i,1,j], N[i,1])
      
      pi[i,1,j] <- lambda[i,1,j]/LAMBDA[i,1]
      
      lambda[i,1,j] <- beta[i,1,j] * z[j]
      
      log(beta[i,1,j]) <- beta0[i,1] + beta1[i] * cov.st[1,j]
      
      # for(k in 1:nreps){
      #   y[i,1,j,k] ~ dbin(p, n[i,1,j])
      # }#end k
    }#end j
      
    for(t in 2:nyears){
      
      beta0[i,t] ~ dnorm(0, 0.01)
      
      log(gamma[i,t-1]) <- gamma0[i] + gamma1[i] * cov.t[t-1]
      
      # for(k in 1:nreps){
      #   Y[i,t,k] ~ dbin(p, N[i,t])
      # }#end k
      
      N[i,t] ~ dpois(LAMBDA[i,t])
      
      LAMBDA[i,t] <- LAMBDA[i,t-1] * gamma[i,t-1]
      
      # n[i,t,1:nsites] ~ dmulti(pi[i,t,1:nsites], N[i,t])
      # 
      # pi[i,t,1:nsites] <- lambda[i,t,1:nsites]/LAMBDA[i,t]
      
      for(j in 1:nsites){
        
        n[i,t,j] ~ dbin(pi[i,t,j], N[i,t])
        
        pi[i,t,j] <- lambda[i,t,j]/LAMBDA[i,t]
        
        lambda[i,t,j] <- beta[i,t,j] * z[j]
        
        log(beta[i,t,j]) <- beta0[i,t] + beta1[i] * cov.st[t,j]
        
        # for(k in 1:nreps){
        #   y[i,t,j,k] ~ dbin(p, n[i,t,j])
        # }

        # y[i,t,j] ~ dpois(p * n[i,t,j])
        
      }#end j
        
    }#end t
      
  }#end i

})
      
#-Compile data-#

data <- list(#y = y,
             #Y = apply(y, MARGIN = c(1,2,4), sum),
             n = n,
             N = N,
             cov.t = cov.t,
             cov.st = cov.st,
             cov.use = cov.use)

constants <- list(nspecies = nspecies, nyears = nyears, nsites = nsites)

#-Initial values-#

# beta0.fun <- function(){
#   
# }
# 
# gamma0.fun <-
#   
# gamma1.fun <-

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
                         z = z
                         )}

#-Parameters to save-#

params <- c(
  "beta0",
  "beta1",
  "gamma0",
  "gamma1",
  #"LAMBDA"
  #"p"
  "omega0", 
  "omega1"
)

params2 <- c("N")

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
nt2 <- 10

#-Run model-#

out1 <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt, #thin2 = nt2,
               samplesAsCodaMCMC = TRUE)
      
