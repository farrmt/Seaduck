#-Library-#
begin <- Sys.time()

library(nimble)
library(coda)
library(parallel)

#-Load data-#
load("./DataFormatting/FormattedData_10km.Rds")
load("./DataFormatting/FormattedCon_10km.Rds")
load("./DataFormatting/FormattedCameraData.Rds")

#-Integrated Community Model-#

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
  
  mu.gamma2 ~ dnorm(0, 0.01)
  tau.gamma2 ~ dgamma(1, 0.001)
  sigma2.gamma2 <- 1/tau.gamma2
  
  tau.gamma.basin ~ dgamma(1, 0.001)
  sigma2.gamma.basin <- 1/tau.gamma.basin
  
  #Population density
  mu.beta0 ~ dnorm(0, 0.01)
  tau.beta0 ~ dgamma(1, 0.001)
  sigma2.beta0 <- 1/tau.beta0
  
  mu.beta1 ~ dnorm(0, 0.01)
  tau.beta1 ~ dgamma(1, 0.001)
  sigma2.beta1 <- 1/tau.beta1
  
  mu.beta2 ~ dnorm(0, 0.01)
  tau.beta2 ~ dgamma(1, 0.001)
  sigma2.beta2 <- 1/tau.beta2
  
  mu.beta3 ~ dnorm(0, 0.01)
  tau.beta3 ~ dgamma(1, 0.001)
  sigma2.beta3 <- 1/tau.beta3
  
  #SPECIES-SPECIFIC PARAMETERS
  
  for(i in 1:nspecies){
    
    for(b in 1:nbasins){
      gamma0[i,b] ~ dnorm(mu.g0[i], tau.gamma.basin)
    }
    
    #Population growth
    mu.g0[i] ~ dnorm(mu.gamma0, tau.gamma0)
    gamma1[i] ~ dnorm(mu.gamma1, tau.gamma1) #PDO
    gamma2[i] ~ dnorm(mu.gamma2, tau.gamma2) #NPGO
    
    #Population density
    beta0[i] ~ dnorm(mu.beta0, tau.beta0) #intercept
    beta1[i] ~ dnorm(mu.beta1, tau.beta1) #depth
    beta2[i] ~ dnorm(mu.beta2, tau.beta2) #depth quadratic
    beta3[i] ~ dnorm(mu.beta3, tau.beta3) #sea surface temperature anomalies
    
    #Correction
    log.correction[i] <- log(correction[i])
    
    #------------#
    #-LIKELIHOOD-#
    #------------#
    
    for(j in 1:nsites[1]){
      #for(j in 1:nsites){
      
      #rho[i,1,j] ~ dgamma(zeta[i], zeta[i])
      rho[i,1,j] ~ dgamma(zeta.bar, zeta.bar)
      
      y[i,1,j] ~ dpois(lambda[i,1,j] * rho[i,1,j])
      
      lambda[i,1,j] <- exp(beta0[i] + 
                             beta1[i] * depth.site[1,j] + 
                             beta2[i] * depth.site[1,j] * depth.site[1,j] +
                             beta3[i] * ssta.site[1,j] +
                             log(area[1,j]) + log.correction[i])
      
    }#end j
    
    for(t in 2:nyears){
      
      # for(b in 1:nbasins){
      # 
      #   log(gamma[i,b,t-1]) <- gamma0[i,b] +
      #     gamma1[i] * PDO[t-1]
      # 
      # }#end b
      
      for(j in 1:nsites[t]){
        
        rho[i,t,j] ~ dgamma(zeta.bar, zeta.bar)
        
        y[i,t,j] ~ dpois(lambda[i,t,j] * rho[i,t,j])
        
        # lambda[i,t,j] <- exp(beta0[i] +
        #                        log(prod(gamma[i,basin[t,j],1:(t-1)])) +
        #                        beta1[i] * depth.site[t,j] +
        #                        log(area[t,j]))
        
        lambda[i,t,j] <- exp(beta0[i] +
                               beta1[i] * depth.site[t,j] +
                               beta2[i] * depth.site[t,j] * depth.site[t,j] +
                               beta3[i] * ssta.site[t,j] +
                               log(area[t,j]) +
                               log.correction[i] +
                               gamma0[i, basin[t,j]] * (year[t] - 1) +
                               gamma1[i] * PDO.c[t-1] +
                               gamma2[i] * NPGO.c[t-1])
      }#end j
      
    }#end t
    
  }#end i
  
  #-Camera data-#
  
  #Composition as proportion of species-specific abundance to community-abundance
  #pi[1:(nspecies-1)] <- theta[1:(nspecies-1)]/THETA

  
  #Linear predictor of species-specific abundance
  for(i in 1:6){
    
    exp.beta0[i] <- exp(beta0[i])
    pi[i] <- exp(beta0[i])/sum(exp.beta0[1:(nspecies-1)])
    
    log(theta[i]) <- beta0[i] + mu.g0.star[i] * 11 + gamma1[i] * PDO.c[11] + gamma2[i] * NPGO.c[11] + log(0.2) 
    pi.star[i] <- theta[i]/THETA
    phi.star[i] <- (phi[i] * sum(exp.beta0[1:(nspecies-1)]) * exp(mu.g0.star[i] * 11 + gamma1[i] * PDO.c[11] + gamma2[i] * NPGO.c[11] + log(0.2)))/THETA
    
    mu.g0.star[i] ~ dnorm(mu.g0[i], tau.gamma.basin)
  }
  for(i in 7:(nspecies-1)){
    
    exp.beta0[i] <- exp(beta0[i+1])
    pi[i] <- exp(beta0[i+1])/sum(exp.beta0[1:(nspecies-1)])
    
    log(theta[i]) <- beta0[i+1] + mu.g0.star[i] * 11 + gamma1[i+1] * PDO.c[11] + gamma2[i+1] * NPGO.c[11] + log(0.2)
    pi.star[i] <- theta[i]/THETA
    phi.star[i] <- (phi[i] * sum(exp.beta0[1:(nspecies-1)]) * exp(mu.g0.star[i] * 11 + gamma1[i+1] * PDO.c[11] + gamma2[i+1] * NPGO.c[11] + log(0.2)))/THETA
    
    mu.g0.star[i] ~ dnorm(mu.g0[i+1], tau.gamma.basin)
  }
  
  for(k in 1:nreps){
    
    rho.bar[k] ~ dgamma(zeta.bar, zeta.bar)
    
    #Front facing camera (true) composition
    FF[1:(nspecies-1),k] ~ dmulti(pi.star[1:(nspecies-1)], FF.total[k])
    
    #Front facing camera (true) community-abundance
    FF.total[k] ~ dpois(THETA * rho.bar[k])
    
    
    for(o in 1:nobs){
      
      #Observer composition
      OBS[1:(nspecies-1),k,o] ~ dmulti(phi.star[1:(nspecies-1)], OBS.total[k,o])
      
      #Observer community-abundance
      OBS.total[k,o] ~ dpois(THETA * E.epsilon[k,o] * rho.bar[k])
      
      #Linear predictor of observation error
      log(E.epsilon[k,o]) <- int.epsilon + alpha * seat[k,o]
      # logit(E.epsilon[k,o]) <- logit(int.epsilon) + alpha * seat[k,o]
      
    }#end o
    
  }#end k
  
  #Correction factor for species-specific observations
  correction[1:6] <- exp(int.epsilon) * phi[1:6]/pi[1:6]
  correction[7] <- exp(int.epsilon)
  correction[8:11] <- exp(int.epsilon) * phi[7:10]/pi[7:10]
  
  phi.ones[1:(nspecies-1)] <- 1
  phi[1:(nspecies-1)] ~ ddirch(phi.ones[1:(nspecies-1)])
  
  #Intercept on observation error
  int.epsilon ~ dnorm(0, 0.01)
  #int.epsilon ~ dunif(0, 1)
  
  #Effect of rear seat
  alpha ~ dnorm(0, 0.01)
  
  #Community-wide expected abundance
  THETA <- sum(theta[1:(nspecies-1)])

  zeta.bar ~ dgamma(4, 2)
  
})

#-Compile data-#

data <- list(y = data.list$n,
             FF = camera.list$FF[-7,], FF.total = camera.list$FF.total,
             OBS = camera.list$OBS[-7,,], OBS.total = camera.list$OBS.total
)

con <- list(nspecies = con.list$nspecies, nyears = con.list$nyears, nsites = con.list$nsites,
            area = con.list$area,
            PDO.c = con.list$PDO.c,
            NPGO.c = con.list$NPGO.c,
            year = con.list$year,
            depth.site = con.list$depth.site,
            ssta.site = con.list$ssta.site,
            nbasins = con.list$nbasins, basin = con.list$basin,
            nreps = camera.list$nreps, nobs = camera.list$nobs,
            seat = camera.list$seat)

#-Initial values-#

beta0.fun <- function(){
  beta0 <- NULL
  beta0[1] <- runif(1, -2.5, -1.5)
  beta0[2] <- runif(1, -3, -1)
  beta0[3] <- runif(1, 0.7, 0.9)
  beta0[4] <- runif(1, -1.2, -0.8)
  beta0[5] <- runif(1, -4, 0)
  beta0[6] <- runif(1, -0.5, 0.5)
  beta0[7] <- runif(1, -2.6, -2.3)
  beta0[8] <- runif(1, -0.6, -0.4)
  beta0[9] <- runif(1, 0.2, 0.6)
  beta0[10] <- runif(1, 1.2, 1.4)
  beta0[11] <- runif(1, 1.3, 1.5)
  return(beta0)
}

beta1.fun <- function(){
  beta1 <- NULL
  beta1[1] <- runif(1, -2.3, -2.1)
  beta1[2] <- runif(1, -2.2, -2)
  beta1[3] <- runif(1, -1.8, -1.7)
  beta1[4] <- runif(1, -1.9, -1.8)
  beta1[5] <- runif(1, -1.7, -1.6)
  beta1[6] <- runif(1, -1.5, -1.4)
  beta1[7] <- runif(1, -2.6, -2.4)
  beta1[8] <- runif(1, -0.7, -0.6)
  beta1[9] <- runif(1, -1.3, -1.2)
  beta1[10] <- runif(1, -1.2, -1)
  beta1[11] <- runif(1, -1.2, -1)
  return(beta1)
}

beta2.fun <- function(){
  beta2 <- NULL
  beta2[1] <- runif(1, 0.1, 0.2)
  beta2[2] <- runif(1, 0.4, 0.5)
  beta2[3] <- runif(1, 0.18, 0.24)
  beta2[4] <- runif(1, 0.2, 0.3)
  beta2[5] <- runif(1, 0.1, 0.2)
  beta2[6] <- runif(1, -0.05, 0.05)
  beta2[7] <- runif(1, -0.8, -0.2)
  beta2[8] <- runif(1, -0.2, -0.1)
  beta2[9] <- runif(1, 0.15, 0.2)
  beta2[10] <- runif(1, 0.06, 0.12)
  beta2[11] <- runif(1, 0.08, 0.14)
  return(beta2)
}


beta3.fun <- function(){
  beta3 <- NULL
  beta3[1] <- runif(1, -0.4, -0.3)
  beta3[2] <- runif(1, 0.35, 0.45)
  beta3[3] <- runif(1, 0.05, 0.15)
  beta3[4] <- runif(1, 0.25, 0.35)
  beta3[5] <- runif(1, 0.25, 0.35)
  beta3[6] <- runif(1, -0.1, 0)
  beta3[7] <- runif(1, -0.35, -0.2)
  beta3[8] <- runif(1, -0.5, -0.4)
  beta3[9] <- runif(1, 0.15, 0.25)
  beta3[10] <- runif(1, 0.18, 0.24)
  beta3[11] <- runif(1, 0.15, 0.25)
  return(beta3)
}

gamma1.fun <- function(){
  gamma1 <- NULL
  gamma1[1] <- runif(1, 0.11, 0.13)
  gamma1[2] <- runif(1, -0.03, 0.01)
  gamma1[3] <- runif(1, 0, 0.02)
  gamma1[4] <- runif(1, 0.22, 0.26)
  gamma1[5] <- runif(1, 0.08, 0.11)
  gamma1[6] <- runif(1, 0.06, 0.1)
  gamma1[7] <- runif(1, 0.06, 0.1)
  gamma1[8] <- runif(1, 0, 0.04)
  gamma1[9] <- runif(1, 0.04, 0.06)
  gamma1[10] <- runif(1, 0, 0.02)
  gamma1[11] <- runif(1, 0.02, 0.05)
  return(gamma1)
}

gamma2.fun <- function(){
  gamma2 <- NULL
  gamma2[1] <- runif(1, 0.12, 0.14)
  gamma2[2] <- runif(1, -0.18, -0.15)
  gamma2[3] <- runif(1, 0.05, 0.06)
  gamma2[4] <- runif(1, 0.14, 0.16)
  gamma2[5] <- runif(1, 0.09, 0.11)
  gamma2[6] <- runif(1, 0.03, 0.05)
  gamma2[7] <- runif(1, 0.02, 0.05)
  gamma2[8] <- runif(1, 0.05, 0.07)
  gamma2[9] <- runif(1, 0.06, 0.08)
  gamma2[10] <- runif(1, 0.02, 0.04)
  gamma2[11] <- runif(1, 0.02, 0.04)
  return(gamma2)
}

pi.init <- apply(camera.list$FF[-7,], 1, sum)/sum(camera.list$FF.total)
#phi.init <- apply(camera.list$OBS[-7,,], 1, sum)/sum(camera.list$OBS.total)

inits <- function(){list(
  mu.gamma0 = runif(1, 0, 0.1),
  sigma2.gamma0 = runif(1, 0, 0.01),
  sigma2.gamma.basin = runif(1, 0, 0.01),
  mu.gamma1 = runif(1, -0.1, 0.1),
  sigma2.gamma1 = runif(1, 0, 0.02),
  mu.gamma2 = runif(1, -0.1, 0.1),
  sigma2.gamma2 = runif(1, 0, 0.02),
  # sigma2.year = runif(1, 0, 0.5),
  mu.beta0 = runif(1, -1, 1),
  sigma2.beta0 = runif(1, 0, 1),
  mu.beta1 = runif(1, -2, -1),
  sigma2.beta1 = runif(1, 0, 0.1),
  mu.beta2 = runif(1, 0, 0.2),
  sigma2.beta2 = runif(1, 0, 0.1),
  mu.beta3 = runif(1, 0, 0.2),
  sigma2.beta3 = runif(1, 0, 0.1),
  # sigma2.beta.basin = runif(1, 1, 3),
  # eps.year = eps.year.fun(),
  # beta.basin = beta.basin.fun(),
  # gamma0 = gamma0.fun(),
  gamma1 = gamma1.fun(),
  gamma2 = gamma2.fun(),
  beta0 = beta0.fun(),
  beta1 = beta1.fun(),
  beta2 = beta2.fun(),
  beta3 = beta3.fun(),
  pi = pi.init,
  int.epsilon = runif(1, -1, 1),
  alpha = runif(1, -0.15, 0.15),
  zeta.bar = runif(1, 0.2, 0.3)
)}

#-Parameters to save-#

params <- c(
  "mu.gamma0",
  "sigma2.gamma0",
  "mu.gamma1",
  "sigma2.gamma1",
  "mu.gamma2",
  "sigma2.gamma2",
  "mu.beta0",
  "sigma2.beta0",
  "mu.beta1",
  "sigma2.beta1",
  "mu.beta2",
  "sigma2.beta2",
  "mu.beta3",
  "sigma2.beta3",
  "sigma2.gamma.basin",
  "int.epsilon",
  "alpha",
  "zeta.bar",
  "beta0",
  "beta1",
  "beta2",
  "beta3",
  "gamma1",
  "gamma2",
  "mu.g0",
  "mu.g0.star",
  "pi",
  "phi",
  "correction")


params2 <- c("gamma0")

#-Parallel function-#

par.fun <- function(seed, data, code, con, params, params2, inits){
  
  library(nimble)
  
  #-MCMC settings-# #10.365 Hours #10 km 5.76 mins, 13.99 mins
  begin <- Sys.time()
  model <- nimbleModel(code = code,
                       constants = con,
                       inits = inits,
                       data = data)
  
  print(Sys.time() - begin)
  
  begin <- Sys.time()
  MCMCconf <- configureMCMC(model, monitors = params, monitors2 = params2)
  print(Sys.time() - begin) #6.78 Days #10 km 1.96 hours
  
  MCMCconf$removeSampler(target = c("int.epsilon", "alpha", "mu.g0.star"))
  
  MCMCconf$addSampler(target = c("int.epsilon", "alpha"),
                      type = "AF_slice")
  
  MCMCconf$addSampler(target = "mu.g0.star",
                      type = "AF_slice")
  
  for(i in 1:6){
    MCMCconf$addSampler(target = c(paste0("phi[", i, "]"), paste0("beta0[", i, "]")),
                        type = "AF_slice")
    MCMCconf$addSampler(target = c(paste0("phi[", i, "]"), paste0("gamma1[", i, "]")),
                        type = "AF_slice")
    MCMCconf$addSampler(target = c(paste0("mu.g0.star[", i, "]"), paste0("beta0[", i, "]")),
                        type = "AF_slice")
  }
  
  for(i in 7:10){
    MCMCconf$addSampler(target = c(paste0("phi[", i, "]"), paste0("beta0[", i+1, "]")),
                        type = "AF_slice")
    MCMCconf$addSampler(target = c(paste0("phi[", i, "]"), paste0("gamma1[", i+1, "]")),
                        type = "AF_slice")
    MCMCconf$addSampler(target = c(paste0("mu.g0.star[", i, "]"), paste0("beta0[", i+1, "]")),
                        type = "AF_slice")
  }
  
  # for(i in 1:10){
  #   MCMCconf$addSampler(target = c(paste0("phi[", i, "]"), paste0("mu.g0.star[", i, "]")),
  #                       type = "AF_slice")
  # }
  
  
  #MCMCconf$removeSampler(target = c("gamma0", "gamma1"))
  # 
  # for(i in 1:(con.list$nspecies - 1)){
  #   MCMCconf$addSampler(target = c(paste0("gamma0[", i, "]"), paste0("gamma1[", i, "]")),
  #                       type = "AF_slice")
  # }
  
  # MCMCconf$removeSampler(target = c("beta0", "mu.b0"))
  # 
  # for(i in 1:(con.list$nspecies - 1)){
  #   MCMCconf$addSampler(target = c(paste0("beta0[", i, ",", 1:9, "]"), paste0("mu.b0[", i, "]")),
  #                       type = "AF_slice")
  # }
  
  begin <- Sys.time()
  MCMC <- buildMCMC(MCMCconf)                                                                                                                                                                                              
  print(Sys.time() - begin) #10 km 35.36 min
  
  begin <- Sys.time()
  compiled.model <- compileNimble(model, MCMC)
  print(Sys.time() - begin) #10 km 56.92 mins
  
  
  
  # nc <-1
  # ni <- 40000
  # nb <- 20000
  # nt <- 20
  
  nc <- 1
  ni <- 13000
  nb <- 3000
  nt <- 10
  
  #-Run model-#
  begin <- Sys.time()
  out.integrated <- runMCMC(compiled.model$MCMC,
                            niter = ni, nburnin = nb,
                            nchains = nc, thin = nt, thin2 = nt,
                            samplesAsCodaMCMC = TRUE)
  print(Sys.time() - begin) #10 km @ 3 chains: 7.33 hours
  #-Save output-#
  return(out.integrated)
}

#-Run model-#

clustID <- makeCluster(3)

out.integrated <- parLapply(cl = clustID, X = 1:3,
                            fun = par.fun,
                            code = code,
                            data = data,
                            con = con,
                            inits = inits(),
                            params = params,
                            params2 = params2)

stopCluster(clustID)

#-Save output-#

out.integrated <- list(samples = mcmc.list(out.integrated[[1]]$samples, out.integrated[[2]]$samples, out.integrated[[3]]$samples), 
                       samples2 = mcmc.list(out.integrated[[1]]$samples2, out.integrated[[2]]$samples2, out.integrated[[3]]$samples2))

# ID <- paste("chain", length(list.files(pattern = "chain", full.names = FALSE)) + 1, sep="")
# assign(ID, out)
# save(list = ID, file = paste0(ID, ".Rds"))
save(out.integrated, file = "./DataAnalysis/integratedmodel_out3.Rds")

end <- Sys.time()

print(end - begin)
