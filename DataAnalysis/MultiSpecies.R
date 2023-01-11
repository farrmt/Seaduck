#-Library-#
begin <- Sys.time()

library(nimble)
library(coda)

#-Load data-#
load("./DataFormatting/FormattedData.Rds")
load("./DataFormatting/FormattedCon.Rds")

#-Nimble code-#

code <- nimbleCode({
  
  #--------#
  #-PRIORS-#
  #--------#
  
  #HYPERPARAMETERS
  
  #Population growth
  mu.gamma0 ~ dnorm(0, 0.01)
  tau.gamma0 ~ dgamma(1, 0.001)
  sigma2.gamma0 <- 1/tau.gamma0
  
  mu.gamma1 ~ dnorm(0, 0.01)
  tau.gamma1 ~ dgamma(1, 0.001)
  sigma2.gamma1 <- 1/tau.gamma1
  
  tau.eps.year ~ dgamma(1, 0.001)
  sigma2.eps.year <- 1/tau.eps.year
  
  #Population density
  mu.beta0 ~ dnorm(0, 0.01)
  tau.beta0 ~ dgamma(1, 0.001)
  sigma2.beta0 <- 1/tau.beta0
  
  mu.beta1 ~ dnorm(0, 0.01)
  tau.beta1 ~ dgamma(1, 0.001)
  sigma2.beta1 <- 1/tau.beta1
  
  tau.eps.basin ~ dgamma(1, 0.001)
  sigma2.eps.basin <- 1/tau.eps.basin
  
  #Zero inflation parameter
  mu.omega0 ~ dunif(0, 1)
  mu.omega0L <- logit(mu.omega0)
  tau.omega0 ~ dgamma(1, 0.001)
  sigma2.omega0 <- 1/tau.omega0
  
  mu.omega1 ~ dnorm(0, 0.01)
  tau.omega1 ~ dgamma(1, 0.001)
  sigma2.omega1 <- 1/tau.omega1

  mu.omega2 ~ dnorm(0, 0.01)
  tau.omega2 ~ dgamma(1, 0.001)
  sigma2.omega2 <- 1/tau.omega2
  
  #RANDOM EFFECTS
  
  for(t in 1:(nyears-1)){
    eps.year[t] ~ dnorm(0, tau.eps.year)
  }
  
  for(b in 1:nbasins){
    eps.basin[b] ~ dnorm(0, tau.eps.basin)
  }
  
  #SPECIES-SPECIFIC PARAMETERS
  
  for(s in 1:nspecies){
    
    #Population growth
    gamma0[s] ~ dnorm(mu.gamma0, tau.gamma0)
    gamma1[s] ~ dnorm(mu.gamma1, tau.gamma1)
    
    #Population density
    beta0[s] ~ dnorm(mu.beta0, tau.beta0) #intercept
    beta1[s] ~ dnorm(mu.beta1, tau.beta1) #depth
    
    #Occupancy
    omega0[s] ~ dnorm(mu.omega0L, tau.omega0) #intercept
    omega1[s] ~ dnorm(mu.omega1, tau.omega1) #auto-regressive
    omega2[s] ~ dnorm(mu.omega2, tau.omega2) #depth
  
  #------------#
  #-LIKELIHOOD-#
  #------------#
    
    N[s,1] ~ dpois(LAMBDA[s,1] * effort[1])
    
    LAMBDA[s,1] ~ dunif(0,10000)
    
    for(j in 1:nsites[1]){
      
      n[s,1,j] ~ dbin(pi[s,1,j], N[s,1])
      
      pi[s,1,j] <- lambda[s,1,j]/(LAMBDA[s,1] * effort[1])
      
      lambda[s,1,j] <- z[s,1,hexagon[1,j]] * exp(beta0[s] +
                                                 beta1[s] * depth.site[1,j] +
                                                 eps.basin[basin[1,j]] +
                                                 log(area[1,j]))
      
    }#end j
      
    for(t in 2:nyears){
      
      log(gamma[s,t-1]) <- gamma0[s] + gamma1[s] * PDO[t-1] + eps.year[t-1]
      
      N[s,t] ~ dpois(LAMBDA[s,t] * effort[t])
      
      LAMBDA[s,t] <- LAMBDA[s,t-1] * gamma[s,t-1]
      
      for(j in 1:nsites[t]){
        
        n[s,t,j] ~ dbin(pi[s,t,j], N[s,t])
        
        pi[s,t,j] <- lambda[s,t,j]/(LAMBDA[s,t] * effort[t])
        
        lambda[s,t,j] <- z[s,t,hexagon[t,j]] * exp(beta0[s] +
                                                   log(prod(gamma[s,1:(t-1)])) +
                                                   beta1[s] * depth.site[1,j] +
                                                   eps.basin[basin[1,j]] +
                                                   log(area[1,j]))
        
      }#end j
      
      for(h in 1:nhex){
        
        z[s,t,h] ~ dbern(omega[s,t,h])
        
        logit(omega[s,t,h]) <- omega0[s] +
                               omega1[s] * z[s,t-1,h] +
                               omega2[s] * depth.hex[h]
        
      }#end h
      
      for(k in 1:npred){
        
        psi[s,t,k] <- 1 - (1 - omega[s,t,hexpred[k]]) ^ COS[k]
        
        pred.lambda[s,t,k] <- exp(beta0[s] + 
                                  log(prod(gamma[s,1:(t-1)])) + 
                                  eps.basin[basinpred[k]] +
                                  beta1[s] * depth.pred[t,k]) * COS[k]
        
        density[s,t,k] <- psi[s,t,k] * pred.lambda[s,t,k]
        
      }#end k
      
    }#end t
    
    for(h in 1:nhex){
      
      z[s,1,h] ~ dbern(omega[s,1,h])
      
      logit(omega[s,1,h]) <- omega0[s] + 
                             omega2[s] * depth.hex[h]
      
    }#end h
    
    for(k in 1:npred){
      
      psi[s,1,k] <- 1 - (1 - omega[s,1,hexpred[k]]) ^ COS[k]
      
      pred.lambda[s,1,k] <- exp(beta0[s] + 
                                eps.basin[basinpred[k]] +
                                beta1[s] * depth.pred[1,k]) * COS[k]
      
      density[s,1,k] <- psi[s,1,k] * pred.lambda[s,1,k]
      
    }#end k
    
    mean.geometric.growth[s] <- exp(gamma0[s])
    
  }#end s
  
})

#-Compile data-#
data <- list(n = data.list$n[,17:25,], N = data.list$N[,17:25])

con <- list(nspecies = con.list$nspecies, nyears = 9, nsites = con.list$nsites[17:25], nhex = con.list$nhex, nbasins = con.list$nbasins, npred = con.list$npred,
            area = con.list$area[17:25,],
            effort = con.list$effort[17:25], hexagon = con.list$hexagon[17:25,],
            PDO = con.list$PDO[c(18:25,27)], 
            depth.site = con.list$depth.site[17:25,], depth.hex = con.list$depth.hex, depth.pred = con.list$depth.pred[17:25,],
            basin = con.list$basin[17:25,], 
            basinhex = con.list$basinhex,
            basinpred = con.list$basinpred, 
            hexpred = con.list$hexpred, 
            COS = con.list$COS)

#-Initial values-#

eps.year.fun <- function(){
  eps.year <- NULL
  eps.year[1] <- runif(1, -1.2, -1)
  eps.year[2] <- runif(1, -0.7, -0.6)
  eps.year[3] <- runif(1, -0.7, -0.6)
  eps.year[4] <- runif(1, -0.55, -0.5)
  eps.year[5] <- runif(1, -0.1, -0.05)
  eps.year[6] <- runif(1, -0.1, -0.05)
  eps.year[7] <- runif(1, -0.1, 0.1)
  eps.year[8] <- runif(1, -0.2, -0.15)
  return(eps.year)
}

eps.basin.fun <- function(){
  eps.basin <- NULL
  eps.basin[1] <- runif(1, -1.7, -1.68)
  eps.basin[2] <- runif(1, -1.5, -1.45)
  eps.basin[3] <- runif(1, -1.7, -1.6)
  eps.basin[4] <- runif(1, -1.65, -1.6)
  eps.basin[5] <- runif(1, -1.66, -1.6)
  eps.basin[6] <- runif(1, -1.6, -1.5)
  eps.basin[7] <- runif(1, -1.7, -1.68)
  eps.basin[8] <- runif(1, -1.56, -1.5)
  eps.basin[9] <- runif(1, -1.6, -1.55)
  return(eps.basin)
}

beta0.fun <- function(){
  beta0 <- NULL
  beta0[1] <- runif(1, 2.3, 2.35)
  beta0[2] <- runif(1, -1.5, -1)
  beta0[3] <- runif(1, 1.1, 1.2)
  beta0[4] <- runif(1, 1, 1.1)
  beta0[5] <- runif(1, 0.68, 0.74)
  beta0[6] <- runif(1, 0.4, 0.5)
  beta0[7] <- runif(1, -0.5, 0)
  beta0[8] <- runif(1, 0.2, 0.3)
  beta0[9] <- runif(1, 0.6, 0.7)
  beta0[10] <- runif(1, 1.4, 1.45)
  beta0[11] <- runif(1, 2.55, 2.6)
  return(beta0)
}

beta1.fun <- function(){
  beta1 <- NULL
  beta1[1] <- runif(1,-0.09, -0.07)
  beta1[2] <- runif(1, -0.1, 0)
  beta1[3] <- runif(1, -0.06, -0.05)
  beta1[4] <- runif(1, -0.08, -0.06)
  beta1[5] <- runif(1, -0.07, -0.05)
  beta1[6] <- runif(1, -0.08, -0.04)
  beta1[7] <- runif(1, -0.06, 0)
  beta1[8] <- runif(1, -0.08, -0.04)
  beta1[9] <- runif(1, -0.03, 0)
  beta1[10] <- runif(1, -0.12, -0.1)
  beta1[11] <- runif(1, -0.11, -0.08)
  return(beta1)
}

gamma0.fun <- function(){
  gamma0 <- NULL
  gamma0[1] <- runif(1, 0.23, 0.3)
  gamma0[2] <- runif(1, 0.5, 0.55)
  gamma0[3] <- runif(1, 0.48, 0.5)
  gamma0[4] <- runif(1, 0.4, 0.42)
  gamma0[5] <- runif(1, 0.43, 0.44)
  gamma0[6] <- runif(1, 0.4, 0.5)
  gamma0[7] <- runif(1, 0.5, 0.6)
  gamma0[8] <- runif(1, 0.42, 0.43)
  gamma0[9] <- runif(1, 0.49, 0.5)
  gamma0[10] <- runif(1, 0.42, 0.43)
  gamma0[11] <- runif(1, 0.15, 0.16)
  return(gamma0)
}

gamma1.fun <- function(){
  gamma1 <- NULL
  gamma1[1] <- runif(1, -0.1, -0.05)
  gamma1[2] <- runif(1, -0.1, 0)
  gamma1[3] <- runif(1, -0.38, -0.35)
  gamma1[4] <- runif(1, -0.35, -0.3)
  gamma1[5] <- runif(1, -0.3, -0.25)
  gamma1[6] <- runif(1, -0.35, -0.3)
  gamma1[7] <- runif(1, -0.2, -0.1)
  gamma1[8] <- runif(1, -0.36, -0.32)
  gamma1[9] <- runif(1, -0.26, -0.24)
  gamma1[10] <- runif(1, -0.3, -0.28)
  gamma1[11] <- runif(1, -0.4, -0.3)
  return(gamma1)
}

omega0.fun <- function(){
  omega0 <- NULL
  omega0[1] <- runif(1, -1.25, -1.1)
  omega0[2] <- runif(1, -2.4, -2.2)
  omega0[3] <- runif(1, -0.4, -0.2)
  omega0[4] <- runif(1, -0.7, -0.5)
  omega0[5] <- runif(1, -1.5, -1.4)
  omega0[6] <- runif(1, -2, -1.9)
  omega0[7] <- runif(1, -2, -1.8)
  omega0[8] <- runif(1, -2.5, -2.3)
  omega0[9] <- runif(1, -0.5, -0.4)
  omega0[10] <- runif(1, -0.3, -0.2)
  omega0[11] <- runif(1, -1.5, -1.3)
  return(omega0)
}

omega1.fun <- function(){
  omega1 <- NULL
  omega1[1] <- runif(1, 2.2, 2.4)
  omega1[2] <- runif(1, 1.5, 2)
  omega1[3] <- runif(1, 3.9, 4.2)
  omega1[4] <- runif(1, 2.7, 2.9)
  omega1[5] <- runif(1, 1.2, 1.4)
  omega1[6] <- runif(1, 3, 3.2)
  omega1[7] <- runif(1, 1.8, 2.2)
  omega1[8] <- runif(1, 3, 3.3)
  omega1[9] <- runif(1, 2, 2.3)
  omega1[10] <- runif(1, 2.3, 2.6)
  omega1[11] <- runif(1, 1.8, 2.1)
  return(omega1)
}

omega2.fun <- function(){
  omega2 <- NULL
  omega2[1] <- runif(1, -0.6, -0.4)
  omega2[2] <- runif(1, -0.6, -0.4)
  omega2[3] <- runif(1, -1.4, -1.2)
  omega2[4] <- runif(1, -1.1, -0.8)
  omega2[5] <- runif(1, -0.5, -0.3)
  omega2[6] <- runif(1, -0.2, 0)
  omega2[7] <- runif(1, -0.6, -0.4)
  omega2[8] <- runif(1, -0.8, -0.5)
  omega2[9] <- runif(1, -0.75, -0.5)
  omega2[10] <- runif(1, -1, -0.75)
  omega2[11] <- runif(1, -0.7, -0.5)
  return(omega2)
}

inits <- function(){list(
  mu.gamma0 = runif(1, 0.35, 0.45),
  sigma2.gamma0 = runif(1, 0, 0.01),
  mu.gamma1 = runif(1, -0.3, -0.2),
  sigma2.gamma1 = runif(1, 0, 0.02),
  sigma2.eps.year = runif(1, 0, 0.5),
  mu.beta0 = runif(1, 0.5, 1.5),
  sigma2.beta0 = runif(1, 0, 1),
  mu.beta1 = runif(1, -0.08, -0.04),
  sigma2.beta1 = runif(1, 0, 0.001),
  sigma2.eps.basin = runif(1, 1, 3),
  mu.omega0 = runif(1, 0.2, 0.3),
  sigma2.omega0 = runif(1, 0, 1),
  mu.omega1 = runif(1, 2, 2.5),
  sigma2.omega1 = runif(1, 0, 1),
  mu.omega2 = runif(1, -0.7, -0.5),
  sigma2.omega2 = runif(1, 0, 0.1),
  eps.year = eps.year.fun(),
  eps.basin = eps.basin.fun(),
  gamma0 = gamma0.fun(),
  gamma1 = gamma1.fun(),
  beta0 = beta0.fun(),
  beta1 = beta1.fun(),
  omega0 = omega0.fun(),
  omega1 = omega1.fun(),
  omega2 = omega2.fun()
)}

#-Parameters to save-#

params <- c(
  "mu.gamma0",
  "sigma2.gamma0",
  "mu.gamma1",
  "sigma2.gamma1",
  "sigma2.eps.year",
  "mu.beta0",
  "sigma2.beta0",
  "mu.beta1",
  "sigma2.beta1",
  "sigma2.eps.basin",
  "mu.omega0",
  "sigma2.omega0",
  "mu.omega1",
  "sigma2.omega1",
  "mu.omega2",
  "sigma2.omega2",
  "eps.year",
  "eps.basin",
  "gamma0",
  "gamma1",
  "beta0",
  "beta1",
  "omega0",
  "omega1",
  "omega2"
)

# params2 <- c(
#   "density",
#   "mean.geometric.growth"
# )

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = con,
                     inits = inits(),
                     data = data)

MCMCconf <- configureMCMC(model, monitors = params)#, monitors2 = params2)

MCMCconf$removeSampler(target = c("omega0", "omega1", "gamma0", "mu.gamma0",
                                  "gamma1", "mu.gamma1", "eps.year", "beta0",
                                  "eps.basin"))

MCMCconf$addSampler(target = c("omega0", "omega1"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("gamma0", "mu.gamma0"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("gamma1", "mu.gamma1"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("gamma0", "eps.year"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("beta0", "eps.basin"),
                    type = "AF_slice")


MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 1
ni <- 20000
nb <- 10000
#nb<- 1000
nt <- 10
# nt2 <- 20

#-Run model-#

out1 <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt, #thin2 = nt2,
               samplesAsCodaMCMC = TRUE)

#-Save output-#

ID <- paste("chain", length(list.files(pattern = "chain", full.names = FALSE)) + 1, sep="")
assign(ID, out)
save(list = ID, file = paste0(ID, ".Rds"))

end <- Sys.time()

print(end - begin)
