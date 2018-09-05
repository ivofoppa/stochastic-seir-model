## Discrete stochastic transmission model with contact matrix
## Latent period (from Comflu_bas): 30% 1 day; 50% 2 days and 20% 3 days;
## Infectious period: 30% 3, 40% 4, 20% 5, 10% 6 days
## This is the POLYMOD matrix for UK all contacts, regrouped
populationData <- read.table('pop_reg5.dat',header = T)

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

latentPeriod0 <- 1.5

infectiousPeriod0 <- 2.5
## Construct population:
populationData <- read.table('pop_reg5.dat',header = T)

## ORder of appearance of flu 
## Seed cases of flu, by HHS
seedInf <- rep(10,10)
seedInf[c(2,6)] <- 500
seedInf[c(3,4)] <- 50

durEpidemic <- 300 ## Number of days in epidemic

## September 12 is day 255 of the year--day one of model; assume Maximum R0 on Jan 15 (day 125 of model); 
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

## Function for creating new infections per "type" - latent stage; all age groups
## A list is returned, the first element being the new infections per age group/type, the reduced 
## susceptibles and the newly infected (latent stage)
nsim <- 10

for (sim in 1:nsim){
  R0minls <- runif(10,1.5,2.2)
  R0maxls <- sapply(R0minls, function(x) runif(1,x,2.4))
  
  corrls <- round(runif(10,64-30,64+30)) ## Random peak R0 between ~ Dec 15 and Feb 15

    generationTime <- runif(1,3.8,4.5)
  latentPeriod <- runif(1,latentPeriod0-.5,latentPeriod0 + .5)
  gamma <- 1 - exp(-1/latentPeriod)
  
  infectiousPeriod <- generationTime - latentPeriod
  delta <- 1 - exp(-1/infectiousPeriod)
  
  R0min <- R0minls[k]
  R0max <- R0maxls[k]
  corr <- corrls[k]
  Itotls <- rep(0,durEpidemic)
  for (hhs in 1:10){
    Einit <-  HHSseedInf[[hhs]]
    Iinit <-  Einit
    population <- as.numeric(populationData[hhs,2:5])
    Sinit <- population - Einit - Iinit
    
    # Has the numbers of individuals by "type" by day since infection
    Rinit <-  Sinit*0
    
    inits <- list(Sinit=Sinit,Einit=Einit,Iinit=Iinit,Rinit=Rinit)
    
    time <- 0
    newCases <- NULL
    S <- inits$Sinit
    E <- inits$Einit
    I <- inits$Iinit
    
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
      
      newCases <- c(newCases,sum(newinf))
      
      Sls <- c(Sls,sum(S))
      Els <- c(Els,sum(E))
      Ils <- c(Ils,sum(I))
      
      time <- time + 1
    }
    Itotls <- Itotls + Ils
  }
  assign(paste0('Itotls',k),Itotls)
}

allItotls <- NULL

for (k in 1:nsim){
allItotls <- c(allItotls,eval(parse(text = paste0('Itotls',k))))  
}

sumItotls <- rep(0,length(Itotls1))
for (k in 1:nsim){
  sumItotls <- sumItotls + eval(parse(text = paste0('Itotls',k)))
}
selind <- which(sumItotls/nsim > 1000)

plot(Itotls1[selind], type = 'l',col = 'red',ylim = c(0,max(allItotls)),xlab = 'Outbreak Day',ylab = 'Number Infectious')

cols <- rainbow(nsim)

for (k in 2:nsim){
  eval(parse(text = paste0('lines(Itotls',k,'[selind],col = cols[k], type = \"l\")')))
}
######################################################################################################
  ######################################################################################################
  