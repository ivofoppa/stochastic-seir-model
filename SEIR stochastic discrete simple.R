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

latentPeriod <- 1.5
gamma <- 1 - exp(-1/latentPeriod)

infectiousPeriod <- 2.5
delta <- 1 - exp(-1/infectiousPeriod)
## Construct population:
Npop <- 310e6
population <- rmultinom(1,Npop,populationFractions)

R0 <- 1.5
beta <- R0 / infectiousPeriod / max(eigen(contactMatrix, 
                                              symmetric = FALSE, only.values = TRUE)$values)

parameters <- list(beta=beta,latentPeriod=latentPeriod,infectiousPeriod=infectiousPeriod,
                   population=population,pmat=pmat)

evec <- eigen(contactMatrix)$vectors[,1]
seedInf <- round(evec/sum(evec)*1000) ## number infected per age group  at beginning
#seedInf <- round(populationFractions*1000) ## number infected per age group  at beginning; compatible with SEIR model

durEpidemic <- 300 ## Number of days in epidemic

# Has the numbers of individuals by "type" by day since infection
Einit <-  seedInf
Sinit <- population - Einit

# Has the numbers of individuals by "type" by day since infection
Iinit <-  Sinit*0
Rinit <-  Sinit*0

inits <- list(Sinit=Sinit,Einit=Einit,Iinit=Iinit,Rinit=Rinit)

## Function to collect Earr_list into 
pSE <- function(I){ ## k is age group, I is matrix of # infectious by latent/infectious time type
  lambda <- beta * contactMatrix %*% (I/population)
  return (as.numeric(1-exp(-lambda)))
}

## Function for creating new infections per "type" - latent stage; all age groups
## A list is returned, the first element being the new infections per age group/type, the reduced 
## susceptibles and the newly infected (latent stage)
nsim <- 10
sim <- 0
while ( sim < nsim ){
  time <- 0
  newCases <- NULL
  Sls <- S <- inits$Sinit
  Els <- E <- inits$Einit
  Ils <- I <- inits$Iinit
  
  while (time < durEpidemic){
    pinf <- pSE(I)
    newlat <- rbinom(4,S,pinf)
    newinf <- rbinom(4,E,gamma)
    newrem <- rbinom(4,I,delta)
    
    S <- S - newlat
    E <- E + newlat - newinf
    I <- I + newinf - newrem
    
    newCases <- c(newCases,sum(newinf))
    
    Sls <- rbind(Sls,S)
    Els <- rbind(Els,E)
    Ils <- rbind(Ils,I)
    
    time <- time + 1
  }
  sim <- sim + 1
  assign(paste0('newCases',sim),newCases)
}

cols <- rainbow(nsim)
plot(newCases1[1:150], type = 'l',col = cols[1])

for (k in 2:nsim){
  eval(parse(text = paste0('lines(newCases',k,'[1:150],col = cols[k], type = \"l\")')))
}
######################################################################################################
######################################################################################################
