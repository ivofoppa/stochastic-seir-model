## Discrete stochastic transmission model with contact matrix
## Latent period (from Comflu_bas): 30% 1 day; 50% 2 days and 20% 3 days;
## Infectious period: 30% 3, 40% 4, 20% 5, 10% 6 days
populationFractions <- c(0.06194508, 0.17698821, 0.59678721, 0.16427950) #UK 2011 data, age group ranges: 0-4, 5-19, 20-64, 65+
populationLabels <- c("Ages 0-4", "Ages 5-19", "Ages 20-64", "Ages 65+")
## This is the POLYMOD matrix for UK all contacts, regrouped
PolyMod_matrix <- matrix(nrow = 4, ncol = 4)
PolyMod_matrix[1,] <- c(1.9200000, 0.5203550, 0.4761399, 0.1272683)
PolyMod_matrix[2,] <- c(1.4867476, 8.7525691, 1.7676622, 0.7706686)
PolyMod_matrix[3,] <- c(4.5871960, 5.9603867, 7.8770903, 3.5312914)
PolyMod_matrix[4,] <- c(0.3375179, 0.7153304, 0.9720697, 1.8867659)

## Matrix entry (row i, column j) denotes the number of potentially infectious contacts a single individual
## from group j has with individuals from group i each day. This will be re-scaled to have spectral radius 1.
contactMatrix <- PolyMod_matrix
specrad <- max(eigen(contactMatrix, 
                     symmetric = FALSE, only.values = TRUE)$values)

## Construct population:
Npop <- 310e6
population <- rmultinom(1,Npop,populationFractions)

evec <- eigen(contactMatrix)$vectors[,1]
seedInf <- round(populationFractions*1000) ## number infected per age group  at beginning
#seedInf <- round(populationFractions*1000) ## number infected per age group  at beginning; compatible with SEIR model

durEpidemic <- 300 ## Number of days in epidemic

# Has the numbers of individuals by "type" by day since infection
Einit <-  seedInf
Iinit <- seedInf

# Has the numbers of individuals by "type" by day since infection
Sinit <-  as.integer(population - Iinit - Einit)

inits <- list(Sinit=Sinit,Einit=Einit,Iinit=Iinit)

## Function to collect Earr_list into 
r0 <- function(t,R0min,R0max,corr){
  y <- (sin((t - corr)/365 * 2 * pi ) + 1) / 2 * (R0max - R0min)
  return(y + R0min)
}

beta <- function(t,R0min,R0max,corr) {
  r0(t,R0min,R0max,corr) / infectiousPeriod / specrad
}

# Has the numbers of individuals by "type" by day since infection
## Function to collect Earr_list into 
pSE <- function(t,I,population,R0min,R0max,corr){ ## k is age group, I is matrix of # infectious by latent/infectious time type
  lambda <- beta(t,R0min,R0max,corr) * contactMatrix %*% as.numeric(I/population)
  return (as.numeric(1-exp(-lambda)))
}

R0minls <- runif(nsim,1.9,2.2)
R0maxls <- sapply(R0minls, function(x) runif(1,x,2.2))

corrls <- round(runif(nsim,34-30,34+30)) ## Random peak R0 between ~ Dec 15 and Feb 15
## Function for creating new infections per "type" - latent stage; all age groups
## A list is returned, the first element being the new infections per age group/type, the reduced 
## susceptibles and the newly infected (latent stage)
nsim <- 10
sim <- 0
while ( sim < nsim ){

  generationTime <- runif(1,3.8,4.5)
  latentPeriod <- runif(1,latentPeriod0-.5,latentPeriod0 + .5)
  gamma <- 1 - exp(-1/latentPeriod)
  
  infectiousPeriod <- generationTime - latentPeriod
  delta <- 1 - exp(-1/infectiousPeriod)
  
  R0min <- R0minls[k]
  R0max <- R0maxls[k]
  corr <- corrls[k]
  
  time <- 0

  S <- as.integer(inits$Sinit)
  E <- as.integer(inits$Einit)*0
  I <- as.integer(inits$Iinit)

  Sls <-   NULL
  Els <-   NULL
  Ils <-   NULL
  
  while (time < durEpidemic){
    pinf <- pSE(time,I,population,R0min,R0max,corr)
    
    newlat <- rbinom(4,S,pinf)
    newinf <- rbinom(4,E,gamma)
    newrem <- rbinom(4,I,delta)
    
    S <- S - newlat
    E <- E + newlat - newinf
    I <- I + newinf - newrem
    
    Sls <- c(Sls,sum(S))
    Els <- c(Els,sum(E))
    Ils <- c(Ils,sum(I))
    
    time <- time + 1
  }
  sim <- sim + 1
  assign(paste0('I',sim,'ls'),Ils)
}

allItotls <- NULL

for (k in 1:nsim){
  allItotls <- c(allItotls,eval(parse(text = paste0('I',k,'ls'))))  
}

sumIls <- rep(0,length(I1ls))
for (k in 1:nsim){
  sumIls <- sumIls + eval(parse(text = paste0('I',k,'ls')))
}
selind <- which(sumIls/nsim > 1000)

plot(I1ls[selind], type = 'l',col = 'red',ylim = c(0,max(allItotls)),xlab = 'Outbreak Day',ylab = 'Number Infectious')

cols <- rainbow(nsim)

for (k in 2:nsim){
  eval(parse(text = paste0('lines(I',k,'ls[selind],col = cols[k], type = \"l\")')))
}
######################################################################################################
######################################################################################################
