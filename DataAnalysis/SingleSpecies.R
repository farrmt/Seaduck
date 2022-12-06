#-Library-#
library(nimble)
library(coda)

#-Load data-#
load("./DataFormatting/FormattedData.Rds")
load("./DataFormatting/FormattedCon.Rds")

#-Nimble code-#

code <- nimbleCode({
  
  #PRIORS
  
  #Population growth
  mu.gamma ~ dnorm(0, 0.01)
  tau.gamma ~ dgamma(1, 0.001)
  sigma2.gamma <- 1/tau.gamma
  
  for(t in 1:(nyears-1)){
    gamma0[t] ~ dnorm(mu.gamma, tau.gamma)
  }
  
  gamma1 ~ dnorm(0, 0.01)

  #Population density
  mu.beta ~ dnorm(0, 0.01)
  tau.beta ~ dgamma(1, 0.001)
  sigma2.beta <- 1/tau.beta

  for(b in 1:nbasins){
    beta0[b] ~ dnorm(mu.beta, tau.beta)
  }
  
  # beta0 ~ dnorm(0, 0.01)
  
  beta1 ~ dnorm(0, 0.01) #Depth
  
  #Zero inflation parameter
  omega0 ~ dunif(0.001, 1)
  omega1 ~ dnorm(0, 0.01)
  omega2 ~ dnorm(0, 0.01)
  # omega3 ~ dnorm(0, 0.01)
  # omega4 ~ dnorm(0, 0.01)
    
  #LIKELIHOOD
    
  N[1] ~ dpois(LAMBDA[1] * effort[1])
    
  LAMBDA[1] ~ dunif(0,1000)
    
  for(j in 1:nsites[1]){
      
    n[1,j] ~ dbin(pi[1,j], N[1])
      
    pi[1,j] <- lambda[1,j]/(LAMBDA[1] * effort[1])
    
    lambda[1,j] <- z[1,hexagon[1,j]] * exp(beta0[basin[1,j]] +
                                           beta1 * depth.site[1,j] +
                                           log(area[1,j]))
      
  }#end j
    
  for(t in 2:nyears){
      
    log(gamma[t-1]) <- gamma0[t-1] + gamma1 * PDO[t-1]
    
    N[t] ~ dpois(LAMBDA[t] * effort[t])
      
    LAMBDA[t] <- LAMBDA[t-1] * gamma[t-1]
      
    for(j in 1:nsites[t]){
        
      n[t,j] ~ dbin(pi[t,j], N[t])
        
      pi[t,j] <- lambda[t,j]/(LAMBDA[t] * effort[t])
      
      lambda[t,j] <- z[t,hexagon[t,j]] * exp(beta0[basin[t,j]] +
                                             log(prod(gamma[1:(t-1)])) +
                                             beta1 * depth.site[t,j] +
                                             log(area[t,j]))
 
    }#end j
    
    for(h in 1:nhex){
      z[t,h] ~ dbern(omega[t,h])
      logit(omega[t,h]) <- logit(omega0) +
                                 omega1 * z[t-1,h] + 
                                 omega2 * depth.hex[h] #+
                                 #omega3 * hex.area[t,h] + #should this be added
                                 #omega4 * pred.lambda[t,h]
      
      # pred.lambda[t,h] <- exp(beta0[basinhex[h]] + 
      #                           log(prod(gamma[1:(t-1)])) + 
      #                           beta1 * depth.hex[h]) * COS.hex[h]
      
    }#end h
    
    for(k in 1:npred){
      
      psi[t,k] <- 1 - (1 - omega[t,hexpred[k]]) ^ COS[k]
      
      pred.lambda[t,k] <- exp(beta0[basinpred[k]] + 
                              log(prod(gamma[1:(t-1)])) + 
                              beta1 * depth.pred[t,k]) * COS[k]

      density[t,k] <- psi[t,k] * pred.lambda[t,k]

    }#end k
      
  }#end t
  
  
  for(h in 1:nhex){

    z[1,h] ~ dbern(omega[1,h])
    logit(omega[1,h]) <- logit(omega0) + 
                         omega2 * depth.hex[h] #+
                         #omega3 * hex.area[1,h] +
                         #omega4 * pred.lambda[1,h]
    
    # pred.lambda[1,h] <- exp(beta0[basinhex[h]] + 
    #                         beta1 * depth.hex[h]) * COS.hex[h]

  }#end h
  
  #DERIVED VALUES
  
  mean.geometric.growth <- exp(mu.gamma)
  
  for(k in 1:npred){
    
    psi[1,k] <- 1 - (1 - omega[1,hexpred[k]]) ^ COS[k]
    
    pred.lambda[1,k] <- exp(beta0[basinpred[k]] + 
                             beta1 * depth.pred[1,k]) * COS[k]
    
    density[1,k] <- psi[1,k] * pred.lambda[1,k]
    
  }#end k
  
  
})

#-Compile data-#
data <- list(n = data.list$n[11,,], N = data.list$N[11,])

con <- list(nyears = con.list$nyears, nsites = con.list$nsites, nhex = con.list$nhex, nbasins = con.list$nbasins, npred = con.list$npred,
            area = con.list$area, #hex.area = con.list$hex.area,
            effort = con.list$effort, hexagon = con.list$hexagon,
            PDO = con.list$PDO[3:28], 
            depth.site = con.list$depth.site, depth.hex = con.list$depth.hex, depth.pred = con.list$depth.pred,
            basin = con.list$basin, 
            basinhex = con.list$basinhex,
            basinpred = con.list$basinpred, 
            hexpred = con.list$hexpred, 
            #COS = rep(100, con.list$npred)
            #COS.hex = con.list$COS.hex)
            COS = con.list$COS)

# data <- list(n = data.list$n[5,,], N = data.list$N[5,])
# 
# con <- list(nyears = con.list$nyears, nsites = con.list$nsites, 
#             area = con.list$area, effort = con.list$effort)

#-Initial values-#


# inits <- function(){list(
#   beta0 = runif(7, -3, 1),
#   gamma = runif(6, 0.8, 1.2)
# )}

gamma0.fun <- function(){
  gamma0 <- NULL
  gamma0[1] <- runif(1, 0.4, 0.5)
  gamma0[2] <- runif(1, -0.6, -0.2)
  gamma0[3] <- runif(1, -1, -0.5)
  gamma0[4] <- runif(1, -0.2, 0.2)
  gamma0[5] <- runif(1, 1.5, 1.8)
  gamma0[6] <- runif(1, -1.5, -1)
  gamma0[7] <- runif(1, -0.4, 0)
  gamma0[8] <- runif(1, 0.1, 0.4)
  gamma0[9] <- runif(1, 1.2, 1.4)
  gamma0[10] <- runif(1, 0, 0.2)
  gamma0[11] <- runif(1, -1.6, -1.4)
  gamma0[12] <- runif(1, 0.75, 1)
  gamma0[13] <- runif(1, 0.75, 1)
  gamma0[14] <- runif(1, -1.4, -1.2)
  gamma0[15] <- runif(1, 0.6, 0.7)
  gamma0[16] <- runif(1, -0.5, -0.3)
  gamma0[17] <- runif(1, -0.3, -0.2)
  gamma0[18] <- runif(1, -0.1, 0.1)
  gamma0[19] <- runif(1, 0.35, 0.45)
  gamma0[20] <- runif(1, 0.25, 0.35)
  gamma0[21] <- runif(1, 0.1, 0.2)
  gamma0[22] <- runif(1, -0.6, -0.4)
  gamma0[23] <- runif(1, 0.2, 0.3)
  gamma0[24] <- runif(1, -0.4, -0.2)
  gamma0[25] <- runif(1, -6, -5)
  gamma0[26] <- runif(1, -1, 1)
  return(gamma0)
}

beta0.fun <- function(){
  beta0 <- NULL
  beta0[1] <- runif(1, -2.4, -2.3)
  beta0[2] <- runif(1, -1.9, -1.8)
  beta0[3] <- runif(1, -2.2, -2.1)
  beta0[4] <- runif(1, -1.9, -1.8)
  beta0[5] <- runif(1, -2.8, -2.7)
  beta0[6] <- runif(1, -2.3, -2.2)
  beta0[7] <- runif(1, -1.9, -1.8)
  beta0[8] <- runif(1, -2.1, -2.0)
  beta0[9] <- runif(1, -2.0, -1.9)
  return(beta0)
}

inits <- function(){list(
  mu.beta = runif(1, -3.5, -2.5),
  beta0 = beta0.fun(),
  # beta0 = runif(1, -3.5, -2.5),
  sigma2.beta = runif(1, 1, 2),
  beta1 = runif(1, -1.9, -1.8),
  gamma0 = gamma0.fun(),
  mu.gamma = runif(1, -0.5, 0.5),
  gamma1 = runif(1, -0.5, 0.5),
  #omega = matrix(runif(1, 0.2, 0.3), ncol = nhex, nrow = nyears),
  omega0 = runif(1, 0.2, 0.3),
  omega1 = runif(1, 0.3, 0.5),
  omega2 = runif(1, 0, 0.2),
  #omega3 = runif(1, -1, 1),
  #omega4 = runif(1, -1, 1),
  sigma2.gamma = runif(1, 1, 3)
)}

#-Parameters to save-#

params <- c(
  "beta0",
  "mu.beta",
  "sigma2.beta",
  "beta1",
  "gamma0",
  "mu.gamma",
  "sigma2.gamma",
  "gamma1",
  "omega0",
  "omega1",
  "omega2"
  #"omega3",
  #"omega4"
)

params2 <- c(
  "pred.lambda",
  "psi",
  "density",
  "mean.geometric.growth"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = con,
                     inits = inits(),
                     data = data)

MCMCconf <- configureMCMC(model, monitors = params, monitors2 = params2)

MCMCconf$removeSampler(target = c("omega0", "omega1", "beta0", "mu.beta", "gamma0", "gamma1"))

MCMCconf$addSampler(target = c("omega0", "omega1"),
                    type = "AF_slice")

# MCMCconf$addSampler(target = c("omega2", "omega4"),
#                     type = "AF_slice")

MCMCconf$addSampler(target = c("beta0", "mu.beta", "gamma0[1]"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("gamma0", "gamma1"),
                    type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 20000
#ni <- 10000
nb <- 10000
#nb <- 5000
nt <- 10
nt2 <- 20

#-Run model-#

out <- runMCMC(compiled.model$MCMC,
                niter = ni, nburnin = nb,
                nchains = nc, thin = nt, thin2 = nt2,
                samplesAsCodaMCMC = TRUE)

 
plot(unlist(out$samples[c(1:nc)][,"omega0"]), unlist(out$samples[c(1:nc)][,"omega1"]))
