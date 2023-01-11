#-----------#
#-Libraries-#
#-----------#

library(coda)
library(sf)
library(tidyverse)

#-Functions-#

year.names <- function(x){
  return(format(as.Date("1994", format = "%Y") + lubridate::years(x), "%Y"))
}

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

#-------------------#
#-Load model output-#
#-------------------#

pattern <- "chain"

#Retrospective analysis
files <- list.files(pattern = pattern, full.names = TRUE)

nc <- length(files)

for(i in 1:nc){
  load(files[i])
  print(i)
}

out <- mcmc.list(purrr::map(mget(ls()[grep(pattern, ls())]), 1))
derived <- mcmc.list(purrr::map(mget(ls()[grep(pattern, ls())]), 2))

#rm(list = ls()[grep(pattern, ls())])

#-Load data-#

#hexagon.data <- st_read("./DataFormatting/hexagon.data.shp")

prediction.data <- st_read("./DataFormatting/prediction.data.shp")

load("./DataFormatting/FormattedCon.Rds")

#------------#
#-Parameters-#
#------------#

params <- attr(out[[1]], "dimnames")[[2]]
params2 <- attr(derived[[1]], "dimnames")[[2]]

#---------#
#-Results-#
#---------#

ni <- nc * length(out[[1]][,1])
ni2 <- nc * length(derived[[1]][,1])

nyears <- con.list$nyears
npred <- con.list$npred
nbasins <- con.list$nbasins
basinpred <- con.list$basinpred
hexpred <- con.list$hexpred

psi <- array(NA, dim = c(nyears,npred,ni2))
lambda <- array(NA, dim = c(nyears,npred,ni2))
density <- array(NA, dim = c(nyears,npred,ni2))
LAMBDA <- array(NA, dim = c(nyears,nbasins+1,ni2))
abundance <- array(NA, dim = c(nyears,nbasins+1,ni2))
muloggrowth <- loggrowth <- array(NA, dim = c(nyears-1,ni))

for(tt in 1:nyears){
  target <- grep(paste("density\\[", tt, ",", sep=""), params2)
  target2 <- grep(paste("pred.lambda\\[", tt, ",", sep=""), params2)
  target3 <- grep(paste("psi\\[", tt, ",", sep=""), params2)
  abundance[tt,1,] <- apply(do.call(rbind, derived[c(1:nc)][,target]), MARGIN = 1, sum)
  LAMBDA[tt,1,] <- apply(do.call(rbind, derived[c(1:nc)][,target2]), MARGIN = 1, sum)
  for(k in 2:(nbasins+1)){
    abundance[tt,k,] <- apply(do.call(rbind, derived[c(1:nc)][,target[basinpred == (k-1)]]), MARGIN = 1, sum)
    LAMBDA[tt,k,] <- apply(do.call(rbind, derived[c(1:nc)][,target2[basinpred == (k-1)]]), MARGIN = 1, sum)
  }
  density[tt,,] <- t(do.call(rbind, derived[c(1:nc)][,target]))
  lambda[tt,,] <- t(do.call(rbind, derived[c(1:nc)][,target2]))
  psi[tt,,] <- t(do.call(rbind, derived[c(1:nc)][,target3]))
  if(tt == nyears) break
    loggrowth[tt,] <- unlist(out[c(1:nc)][,grep(paste("gamma0\\[", tt, "\\]", sep = ""), params)]) +
                      unlist(out[c(1:nc)][,grep(paste("gamma1", sep = ""), params)]) * con.list$PDO[tt+2]
    muloggrowth[tt,] <- unlist(out[c(1:nc)][,grep(paste("mu.gamma", sep = ""), params)]) +
                        unlist(out[c(1:nc)][,grep(paste("gamma1", sep = ""), params)]) * con.list$PDO[tt+2]
}

meanLAMBDA <- apply(LAMBDA, MARGIN = c(1,2), mean, na.rm = T)
LAMBDA97.5 <- apply(LAMBDA, MARGIN = c(1,2), quantile, probs = 0.975, na.rm = T)
LAMBDA2.5 <- apply(LAMBDA, MARGIN = c(1,2), quantile, probs = 0.025, na.rm = T)

meanN <- apply(abundance, MARGIN = c(1,2), mean, na.rm = T)
N97.5 <- apply(abundance, MARGIN = c(1,2), quantile, probs = 0.975, na.rm = T)
N2.5 <- apply(abundance, MARGIN = c(1,2), quantile, probs = 0.025, na.rm = T)

meanGrowth <- apply(loggrowth, MARGIN = 1, FUN = function(x){exp(mean(x))})
Growth97.5 <- apply(loggrowth, MARGIN = 1, FUN = function(x){exp(quantile(x, probs = 0.975))})
Growth2.5 <- apply(loggrowth, MARGIN = 1, FUN = function(x){exp(quantile(x, probs = 0.025))})

meanmuGrowth <- apply(muloggrowth, MARGIN = 1, FUN = function(x){exp(mean(x))})
muGrowth97.5 <- apply(muloggrowth, MARGIN = 1, FUN = function(x){exp(quantile(x, probs = 0.975))})
muGrowth2.5 <- apply(muloggrowth, MARGIN = 1, FUN = function(x){exp(quantile(x, probs = 0.025))})

Ndata <- reshape2::melt(meanN, varnames = c("year", "basin"), value.name = "abundance")
Ndata$basin <- factor(Ndata$basin)
#data$year <- as.numeric(year.names(data$year)) #MTF: need to find which years to omit
Ndata$`97.5%` <- reshape2::melt(N97.5, varnames = c("year", "basin"), value.name = "97.5%")$`97.5%`
Ndata$`2.5%` <- reshape2::melt(N2.5, varnames = c("year", "basin"), value.name = "2.5%")$`2.5%`

Gdata <- data.frame(growth = c(meanGrowth, mean(unlist(derived[c(1:nc)][,"mean.geometric.growth"]))),
                    upper = c(Growth97.5, quantile(unlist(derived[c(1:nc)][,"mean.geometric.growth"]), probs = 0.975)),
                    lower = c(Growth2.5, quantile(unlist(derived[c(1:nc)][,"mean.geometric.growth"]), probs = 0.025)),
                    mu.growth = c(meanmuGrowth, mean(unlist(derived[c(1:nc)][,"mean.geometric.growth"]))),
                    mu.upper = c(muGrowth97.5, quantile(unlist(derived[c(1:nc)][,"mean.geometric.growth"]), probs = 0.975)),
                    mu.lower = c(muGrowth2.5, quantile(unlist(derived[c(1:nc)][,"mean.geometric.growth"]), probs = 0.025)),
                    year = (1:nyears),
                    type = c(rep("annual", (nyears-1)), "mean"))

meanD <- apply(density, MARGIN = c(1,2), mean, na.rm = T)#/matrix(rep(con.list$COS, nyears), ncol = npred) #* (10/mean(area, na.rm = T))
meanD <- reshape2::melt(meanD, varnames = c("year", "hexagon"), value.name = "density")
meanD <- data.frame(meanD) %>% arrange(year, hexagon)
Ddata <- do.call(rbind, replicate(nyears, prediction.data, simplify = F))
Ddata$year <- meanD$year
Ddata$density <- meanD$density

meanL <- apply(lambda, MARGIN = c(1,2), mean, na.rm = T)
meanL <- reshape2::melt(meanL, varnames = c("year", "hexagon"), value.name = "lambda")
meanL <- data.frame(meanL) %>% arrange(year, hexagon)
Ldata <- do.call(rbind, replicate(nyears, prediction.data, simplify = F))
Ldata$year <- meanL$year
Ldata$lambda <- meanL$lambda

meanPsi <- apply(psi, MARGIN = c(1,2), mean, na.rm = T)
meanPsi <- reshape2::melt(meanPsi, varnames = c("year", "hexagon"), value.name = "psi")
meanPsi <- data.frame(meanPsi) %>% arrange(year, hexagon)
Psidata <- do.call(rbind, replicate(nyears, prediction.data, simplify = F))
Psidata$year <- meanPsi$year
Psidata$psi <- meanPsi$psi

Figure1A <- ggplot(data = Ndata %>% filter(basin == 1), aes(x = year, y = abundance)) +
  geom_point(col = "grey") +
  geom_line(col = "grey") + 
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.5) +
  #facet_wrap(. ~ basin, scales = 'free', ncol = 2, dir = "v") +
  #scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "italic")) +
  labs(y = "Relative abundance", x = "Year")

Figure1B <- ggplot(data = Gdata %>% filter(type == "annual"), aes(x = year, y = growth)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = Gdata %>% filter(type == "mean") %>% select(mu.growth) %>% .$growth, col = "blue") +
  geom_rect(data = Gdata %>% filter(type == "mean"), aes(ymin = lower, ymax = upper, xmin = 1, xmax = nyears), fill = "blue", alpha = 0.5) +
  geom_point(aes(x = year, y = mu.growth),col = "black", position = position_nudge(x = 0.5)) +
  geom_errorbar(aes(x = year, ymin = mu.lower, ymax = mu.upper), col = "black", width = 0, position = position_nudge(x = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  ylim(0,2) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "italic")) +
  labs(y = "Expected Growth Rate", x = "Year")

Figure1C <- ggplot(data = Gdata %>% filter(type == "annual"), aes(x = year, y = growth)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = Gdata %>% filter(type == "mean") %>% select(growth) %>% .$growth, col = "blue") +
  geom_rect(data = Gdata %>% filter(type == "mean"), aes(ymin = lower, ymax = upper, xmin = 1, xmax = nyears), fill = "blue", alpha = 0.5) +
  geom_point(aes(x = year, y = growth),col = "black", position = position_nudge(x = 0.5)) +
  geom_errorbar(aes(x = year, ymin = lower, ymax = upper), col = "black", width = 0, position = position_nudge(x = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  #xlim(2,27) +
  ylim(0,7) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "italic")) +
  labs(y = "Growth Rate", x = "Year")

cowplot::plot_grid(Figure1A, Figure1B, Figure1C, ncol = 1, align = "v")

Figure3 <- ggplot() +
  geom_sf(data = prediction.data) +
  geom_sf(data = Ddata, aes(fill = density), size = 0.01) + 
  scale_fill_gradient(name = "Abundance", low = "white", high = "black") +
  #facet_wrap(. ~ year) +
  theme_bw() +
  transition_states(year) +
  labs(subtitle = "Year {next_state}")
animate(Figure3)

Figure3 + exit_shrink()

anim_save(file = "PostAnalysis/map.gif")

ggplot() +
  geom_sf(data = prediction.data) +
  geom_sf(data = Ddata, aes(fill = density), size = 0.01) + 
  scale_fill_gradient(low = "white", high = "black") +
  facet_wrap(. ~ year) +
  theme_bw()

ggplot() +
  geom_sf(data = prediction.data) +
  geom_sf(data = Ldata %>% filter(year >= 22), aes(fill = lambda), size = 0.01) + 
  scale_fill_gradient(low = "white", high = "black") +
  facet_wrap(. ~ year) +
  theme_bw()



ggplot() +
  geom_sf(data = prediction.data) +
  geom_sf(data = Psidata %>% filter(year <= 5), aes(fill = psi), size = 0.01) + 
  scale_fill_gradient(low = "green", high = "red") +
  facet_wrap(. ~ year) +
  theme_bw()

#TWS
bath <- as(bathymetry, "SpatialPixelsDataFrame")
bath <- as.data.frame(bath)

ggplot() +
  geom_sf(data = prediction.data, aes(fill = Basin_1), alpha = 0.4, size = 0.001) +
  theme_bw() +
  scale_fill_discrete(name = "Basin")

ggplot() +
  geom_tile(data = bath, aes(x = x, y = y, fill = Bathymetry)) +
  theme_bw() +
  scale_fill_continuous(name = "Depth")

ggplot() +
  geom_sf(data = prediction.data, aes(fill = use), size = 0.1) +
  theme_bw() +
  scale_fill_manual(name = "Presence/Absence", values = c("White", "Grey"))


beta0 <- summary(out)[[1]][1:9,"Mean"]
beta1 <- summary(out)[[1]][10,"Mean"]
gamma0 <- summary(out)[[1]][11:36,"Mean"]
gamma1 <- summary(out)[[1]][37,"Mean"]
omega0 <- summary(out)[[1]]["omega0","Mean"]
omega1 <- summary(out)[[1]]["omega1","Mean"]
omega2 <- summary(out)[[1]]["omega2","Mean"]


gamma <- exp(gamma0 + gamma1 * con.list$PDO[3:28])

loglambda <- matrix(NA, nrow = 7, ncol = 27)
loglambda[1,] <- beta0

prediction.data$density <- NULL

for(i in 1:dim(prediction.data)[1]){
  omega <- expit(logit(omega0) + omega1)
  prediction.data$density[i] <- 
}
density[t,k] <- omega[t,hexpred[k]] * exp(beta0[basinpred[k]] + 
                                            log(prod(gamma[1:(t-1)])) + 
                                            beta1 * depth.pred[t,k]) * COS[k]