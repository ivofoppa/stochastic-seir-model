## Discrete stochastic transmission model with contact matrix
## Latent period (from Comflu_bas): 30% 1 day; 50% 2 days and 20% 3 days;
## Infectious period: 30% 3, 40% 4, 20% 5, 10% 6 days
populationData <- read.table('pop_reg5.dat',header = T)

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

evec <- eigen(contactMatrix)$vectors[,1]
evec <- evec/sum(evec)
## Construct population:
populationData <- read.table('pop_reg5.dat',header = T)
pop <- as.integer(colSums(populationData[,2:5]))
populationFractions <- pop/sum(pop)

## Creating population according to latent (rows) and infectious periods (columns) in each of the age groups
r0 <- function(t,R0min,R0max,corr){
  y <- (sin((t - corr)/365 * 2 * pi ) + 1) / 2 * (R0max - R0min)
  return(y + R0min)
}

beta <- function(t,R0min,R0max,corr,infectiousPeriod) {
  r0(t,R0min,R0max,corr) / infectiousPeriod / specrad
}

pSE <- function(t,I,R0min,R0max,corr,infectiousPeriod){ ## k is age group, I is matrix of # infectious by latent/infectious time type
  lambda <- beta(t,R0min,R0max,corr,infectiousPeriod) * contactMatrix %*% (I/pop)
  return (as.numeric(1-exp(-lambda)))
}

seedInf <- round(evec*100) ## number infected per age group  at beginning
#seedInf <- round(populationFractions*1000) ## number infected per age group  at beginning; compatible with SEIR model

durEpidemic <- 300 ## Number of days in epidemic

# Has the numbers of individuals by "type" by day since infection
Einit <-  seedInf
Iinit <- rep(0,4)

# Has the numbers of individuals by "type" by day since infection
Sinit <-  as.integer(pop - Iinit - Einit)

inits <- list(Sinit=Sinit,Einit=Einit,Iinit=Iinit)

## Function to collect Earr_list into 

## Function for creating new infections per "type" - latent stage; all age groups
## A list is returned, the first element being the new infections per age group/type, the reduced 
## susceptibles and the newly infected (latent stage)
nsim <- 10

simSerialIntervals <- runif(nsim,2.4,2.8)
simLatentPeriods <- runif(nsim,1.3,1.7)
simInfectiousPeriods <- 2*(simSerialIntervals - simLatentPeriods)

R0minls <- runif(nsim,2,2.2)
R0maxls <- sapply(R0minls, function(x) runif(1,x,2.4))

corrls <- round(runif(nsim,64-30,64+30)) ## Random peak R0 between ~ Dec 15 and Feb 15

sim <- 1

while ( sim <= nsim ){
  R0min <- R0minls[sim]
  R0max <- R0maxls[sim]
  corr <- corrls[sim]
  latentPeriod <- simLatentPeriods[sim]
  gamma <- 1/latentPeriod
  infectiousPeriod <- simInfectiousPeriods[sim]
  delta <- 1/infectiousPeriod
  
  
  
  S <- as.integer(inits$Sinit)
  E <- as.integer(inits$Einit)
  I <- as.integer(inits$Iinit)
  
  Sls <-   NULL
  Els <-   NULL
  Ils <-   NULL
  
  time <- 1
  while (time <= durEpidemic){
    pinf <- pSE(time,I,R0min,R0max,corr,infectiousPeriod)
    
    newlat <- sapply(1:4,function(ag) ifelse(S[ag] > 0, rbinom(1,S[ag],pinf[ag]),0))
    newinf <- sapply(1:4,function(ag) ifelse(E[ag] > 0, rbinom(1,E[ag],gamma),0))
    newrem <- sapply(1:4,function(ag) ifelse(I[ag] > 0, rbinom(1,I[ag],delta),0))
    
    S <- S - newlat
    E <- E + newlat - newinf
    I <- I + newinf - newrem
    
    Sls <- c(Sls,sum(S))
    Els <- c(Els,sum(E))
    Ils <- c(Ils,sum(I))
    
    time <- time + 1
  }
  assign(paste0('Infectious',sim),Ils)
  sim <- sim + 1
}

cols <- rainbow(nsim)
xlim <- 110
xlow <- 30
plot(Infectious1[xlow:xlim], type = 'l',col = cols[1],xlab = 'Outbreak Day',ylab = 'Number Infectious')

#plot(Ils[1:150], type = 'l',col = cols[1])

for (k in 2:nsim){
  eval(parse(text = paste0('lines(infectious',k,'[xlow:xlim],col = cols[k], type = \"l\")')))
}
######################################################################################################
######################################################################################################
