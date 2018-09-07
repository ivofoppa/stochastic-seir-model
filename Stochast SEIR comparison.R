## Discrete stochastic transmission model with contact matrix
## Latent period (from Comflu_bas): 30% 1 day; 50% 2 days and 20% 3 days;
## Infectious period: 30% 3, 40% 4, 20% 5, 10% 6 days
populationData <- read.table('pop_HHS_4.dat',header = F)

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
pop <- as.integer(colSums(populationData[,2:5]))
populationFractions <- pop/sum(pop)

## Creating population according to latent (rows) and infectious periods (columns) in each of the age groups

betafn <- function(t,R0,infectiousPeriod) {
  R0 / infectiousPeriod / specrad 
}

pSE <- function(t,I,R0,infectiousPeriod,deltat){ ## k is age group, I is matrix of # infectious by latent/infectious time type
  lambda <- betafn(t,R0,infectiousPeriod) * 1/pop * (contactMatrix %*% I)
  return (as.numeric(1-exp(-lambda*deltat)))
}

seedInfections <- 1000
seedInf <- seedInfections * populationFractions ## number infected per age group  at beginning
#seedInf <- round(populationFractions*1000) ## number infected per age group  at beginning; compatible with SEIR model

durEpidemic <- 300 ## Number of days in epidemic

R0 <- 1.6
latentPeriod <- 1.5
gamma <- 1/latentPeriod
infectiousPeriod <- 2.5 
delta <- 1/infectiousPeriod

# Has the numbers of individuals by "type" by day since infection
Sinit <- round(pop - seedInf)
Einit <- round(pop*0 + seedInf / (1 + gamma / delta))
Iinit <- round(pop*0 + seedInf / (1 + delta / gamma))
Rinit <- round(pop*0 )

inits <- list(Sinit=Sinit,Einit=Einit,Iinit=Iinit,Rinit=Rinit)

## Function to collect Earr_list into 

## Function for creating new infections per "type" - latent stage; all age groups
## A list is returned, the first element being the new infections per age group/type, the reduced 
## susceptibles and the newly infected (latent stage)
Iarr <- newInfarr <- newRemarr <- Sarr <- NULL
deltat <- .1
for (sim in 1:nsim){
  S <- as.integer(inits$Sinit)
  E <- as.integer(inits$Einit)
  I <- as.integer(inits$Iinit)
  R <- as.integer(inits$Rinit)
  
  Sls <- Els <- Ils <- Rls <- newInfls <- newRemls <- NULL
  
  time <- deltat
  while (time <= durEpidemic){
    pinf <- pSE(time,I,R0,infectiousPeriod,deltat)
    
    newlat <- sapply(1:4,function(ag) ifelse(S[ag] > 0, rbinom(1,S[ag],pinf[ag]),0))
    newinf <- sapply(1:4,function(ag) ifelse(E[ag] > 0, rbinom(1,E[ag],1 - exp(- gamma * deltat)),0))
    newrem <- sapply(1:4,function(ag) ifelse(I[ag] > 0, rbinom(1,I[ag],1 - exp(- delta * deltat)),0))
    
    S <- S - newlat
    E <- E + newlat - newinf
    I <- I + newinf - newrem
    R <- R + newrem
    
    Sls <- c(Sls,sum(S))
    Els <- c(Els,sum(E))
    Ils <- c(Ils,sum(I))
    Rls <- c(Rls,sum(R))
    newInfls <- c(newInfls,sum(newinf))
    newRemls <- c(newRemls,sum(newrem))
    
    time <- time + deltat
  }
  Iarr <- rbind(Iarr, Ils,deparse.level = 0)
  Sarr <- rbind(Sarr, Sls,deparse.level = 0)
  newInfarr <- rbind(newInfarr, newInfls,deparse.level = 0)
  newRemarr <- rbind(newRemarr, newRemls,deparse.level = 0)
}

model <- SEIRModel(R0 = R0,
                   latentPeriod = 1.5,
                   infectiousPeriod = 2.5,
                   seedInfections = 1000,
                   simulationLength = 300,priorImmunity = 0,
                   population = Ntot,populationFractions = populationFractions,
                   contactMatrix = contactMatrix)


modelTimes <- model$rawOutput[, "time"]
modelInfections <- getInfectionTimeSeries(model, byGroup = FALSE)
modelSusceptibles <- getSus
selind <- which(modelTimes <= 150)
plot(modelTimes[selind],modelInfections[selind],type = 'l')


#########################################################################################
beta <- betafn(1,R0,infectiousPeriod)
parameters <- c(gamma=gamma,delta=delta,beta=beta)
### Define model
SEIRmod <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    forceOfInfection <- beta / pop * (contactMatrix %*% c(I1,I2,I3,I4))
    
    S_to_E1 <- S1 * forceOfInfection[1]
    E_to_I1 <- gamma * E1
    I_to_R1 <- delta * I1
    
    S_to_E2 <- S2 * forceOfInfection[2]
    E_to_I2 <- gamma * E2
    I_to_R2 <- delta * I2
    
    S_to_E3 <- S3 * forceOfInfection[3]
    E_to_I3 <- gamma * E3
    I_to_R3 <- delta * I3
    
    S_to_E4 <- S4 * forceOfInfection[4]
    E_to_I4 <- gamma * E4
    I_to_R4 <- delta * I4
    
    #Derivatives
    dS1 <- -S_to_E1
    dE1 <- S_to_E1 - E_to_I1
    dI1 <- E_to_I1 - I_to_R1
    dR1 <- I_to_R1
    
    dS2 <- -S_to_E2
    dE2 <- S_to_E2 - E_to_I2
    dI2 <- E_to_I2 - I_to_R2
    dR2 <- I_to_R2
    
    dS3 <- -S_to_E3
    dE3 <- S_to_E3 - E_to_I3
    dI3 <- E_to_I3 - I_to_R3
    dR3 <- I_to_R3
    
    dS4 <- -S_to_E4
    dE4 <- S_to_E4 - E_to_I4
    dI4 <- E_to_I4 - I_to_R4
    dR4 <- I_to_R4
    
    # return the rate of change
    list(c(dS1,dE1,dI1,dR1,
           dS2,dE2,dI2,dR2,
           dS3,dE3,dI3,dR3,
           dS4,dE4,dI4,dR4))
  })
}
### Initial conditions
S10 <- (as.integer(inits$Sinit))[1]
E10 <- (as.integer(inits$Einit))[1]
I10 <- (as.integer(inits$Iinit))[1]
R10 <- (as.integer(inits$Rinit))[1]

S20 <- (as.integer(inits$Sinit))[2]
E20 <- (as.integer(inits$Einit))[2]
I20 <- (as.integer(inits$Iinit))[2]
R20 <- (as.integer(inits$Rinit))[2]

S30 <- (as.integer(inits$Sinit))[3]
E30 <- (as.integer(inits$Einit))[3]
I30 <- (as.integer(inits$Iinit))[3]
R30 <- (as.integer(inits$Rinit))[3]

S40 <- (as.integer(inits$Sinit))[4]
E40 <- (as.integer(inits$Einit))[4]
I40 <- (as.integer(inits$Iinit))[4]
R40 <- (as.integer(inits$Rinit))[4]
###################################################################################################
library(deSolve)
state <- c(S1=S10,E1=E10,I1=I10,R1=R10,
           S2=S20,E2=E20,I2=I20,R2=R20,
           S3=S30,E3=E30,I3=I30,R3=R30,
           S4=S40,E4=E40,I4=I40,R4=R40)
times <- seq(0,durEpidemic,.01)
out <- data.frame(ode(y = state, times = times, func = SEIRmod, parms = parameters))

SEIRIls <- (out$I1 + out$I2 + out$I3 + out$I4)
SEIRSls <- (out$S1 + out$S2 + out$S3 + out$S4)
###################################################################################################

toTime <- 150
selind <- which(times <= toTime)

ymax <- max(c(Iarr,SEIRIls))
# ymax <- max(c(Sarr,SEIRSls))
colls <- rainbow(nsim)

plot(times[selind],SEIRIls[selind],xlab = 'Outbreak Day',ylab = 'Number Infectious',type = 'l', col = 'blue',xlim = c(0,toTime),ylim = c(0,ymax))
for (sim in 1:10){
 lines(seq(deltat,toTime,deltat),Iarr[sim,seq(deltat,toTime,deltat)/deltat],lwd = 1, col = colls[sim])
}

plot(Ils,col = 'green')
#########################################################################################
#########################################################################################
#########################################################################################
library(flumodels)

Ntot <- sum(pop)

model <- SEIRModel(R0 = 1.6,
                   latentPeriod = 1.5,
                   infectiousPeriod = 2.5,
                   seedInfections = seedInfections,
                   simulationLength = 300,priorImmunity = 0,
                   population = Ntot,populationFractions = populationFractions,contactMatrix = contactMatrix)


modelTimes <- model$rawOutput[, "time"]
modelInfections <- getInfectionTimeSeries(model, byGroup = FALSE)