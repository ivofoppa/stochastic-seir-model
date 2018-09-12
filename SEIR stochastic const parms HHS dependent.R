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
## Creating population according to latent (rows) and infectious periods (columns) in each of the age groups
r0 <- function(t,R0){
  return(R0)
}

beta <- function(t,R0,infectiousPeriodMean) {
  r0(t,R0) / infectiousPeriodMean / specrad
}

pSE <- function(t,hhs,HHSIarr_list,HHStypePopulation,R0,infectiousPeriodMean,HHSfact){ ## k is age group, I is matrix of # infectious by latent/infectious time type
  HHSI <- sapply(1:4, function(ag) sum(sapply(HHSIarr_list[-hhs], function(hhsel) sum(unlist(hhsel[[ag]])))))
  I <- sapply(1:4, function(x) sum(unlist(HHSIarr_list[[hhs]][[x]])))
  pop <- sapply(1:4, function(x) sum(HHStypePopulation[[hhs]][[x]]))
  HHSpop <- sapply(1:4, function(ag) sum(sapply(HHStypePopulation[-hhs], function(pop) sum(pop[[ag]]))))
  beta <- beta(t,R0,infectiousPeriodMean)
  lambda <- beta * contactMatrix %*% as.numeric(I/pop) + beta*HHSfact * contactMatrix %*% as.numeric(HHSI/HHSpop) 
  return (as.numeric(1-exp(-lambda)))
}

newInfect <- function(t,hhs,HHSIarr_list,HHSSarr,HHSEarr_list,HHStypePopulation,R0,infectiousPeriodMean,HHSfact){
  newinfList <- list()
  Sarr2 <- HHSSarr[[hhs]]
  Earr_list2 <- HHSEarr_list[[hhs]]
  pinfls <- pSE(t,hhs,HHSIarr_list,HHStypePopulation,R0,infectiousPeriodMean,HHSfact)
  for (ag in 1:4){
    Sarr_ag <- Sarr2[[ag]]
    Earr <- Earr_list2[[ag]]
    pinf <- pinfls[[ag]]
    newinf <- structure(vapply(Sarr_ag, function(x) rbinom(1,x,pinf), numeric(1)), dim=dim(Sarr_ag))
    Earr_list2[[ag]] <- sapply(seq_along(Earr), function(x) {Earr[[x]][,1] <- Earr[[x]][,1] + newinf[,x]; return(Earr[[x]])})
    Sarr2[[ag]] <- Sarr2[[ag]] - newinf
    newinfList[[ag]] <- newinf
  }
  return(list(Sarr2,Earr_list2,newinfList))
}

dayProgI <- function(hhs,HHSEarr_list,HHSIarr_list){
  Iarr_list2 <- HHSIarr_list[[hhs]]
  for (ag in 1:4){
    for (i in 1:4){
      newinf <- sapply(1:3,function(l) HHSEarr_list[[hhs]][[ag]][[l]][i,l])
      Iarr_list2[[ag]][[i]] <- cbind(newinf,Iarr_list2[[ag]][[i]][,1:i-1],deparse.level = 0)
    }
  }
  return(Iarr_list2)
} 

dayProgE <- function(hhs,HHSEarr_list){
  Earr_list2 <- HHSEarr_list[[hhs]]
  for(ag in 1:4){
    for (l in 1:3){
      if (l > 1){
        Earr_list2[[ag]][[l]] <- cbind(rep(0,length(infectiousPeriods)),Earr_list2[[ag]][[l]][,1:(l - 1)])
      } else {
        Earr_list2[[ag]][[l]] <- cbind(rep(0,length(infectiousPeriods)))
      }
    }
  }
  return(Earr_list2)
}


HHStypePopulation <- list()

for (hhs in 1:10){
  ls <- sapply(1:4, function(nag) list(sapply(rmultinom(1,populationData[hhs,nag+1],latentPeriods),function(x) sapply(x, function(y) rmultinom(1,y,infectiousPeriods)))))
  HHStypePopulation[[hhs]] <- ls
}

HHSpmat <- list()
for (hhs in 1:10){
  typePopulation <- HHStypePopulation[[hhs]]
  HHSpmat[[hhs]] <- sapply(1:4,function(x) list(typePopulation[[x]]/sum(typePopulation[[x]])))
}

## Seed cases of flu, by HHS
## Seed cases of flu, by HHS
seedInf <- rep(1,10)
seedInf[c(2,6)] <- 10
seedInf[c(3,4)] <- 2

durEpidemic <- 300
nsim <- 10
R0 <- 2.3
colls <- rainbow(nsim)

latentPeriodMean <- 1.5
SerialInterval <- 2.6  
lp3 <- -1; lp2 <- -1
while (lp3 < 0 | lp3 > 1 | lp2 < 0 | lp2 > 1 ){
  lp1 <- runif(1,.1,.9)
  lp2 <- (latentPeriodMean - 3 + 2*lp1)/(-1)
  lp3 <- 1 - lp1 - lp2
}

latentPeriods <- c(lp1,lp2,lp3) ## frequencies of latent period durations 1-3 days

## Combinations of infectious/latent periods (discrete) are considered "types"
infectiousPeriodMean <- 2*(SerialInterval - latentPeriodMean)

ip4 <- -1; ip3 <- -1
while (ip4 < 0 | ip4 > 1 | ip3 < 0 | ip3 > 1 ){
  ip1 <- runif(1,0,.9)
  ip2 <- runif(1,0,1 - ip1)
  ip3 <- 4 - infectiousPeriodMean - 3*ip1 - 2*ip2
  ip4 <- 1 - ip1 - ip2 - ip3
}

infectiousPeriods <- c(ip1,ip2,ip3,ip4) ## frequencies of infectious durations 3-6 days

HHSfact <- 0.001

univParameters <- list(simSerialIntervals=simSerialIntervals,simLatentPeriods=simLatentPeriods,siminfectiousPeriods=siminfectiousPeriods,HHSfact=HHSfact)
sink('univParameters3.txt')
print(univParameters)
sink()
# Has the numbers of individuals by "type" by day since infection
Itotarr <- ARtotarr <- inctotarr <- NULL

HHSEarr_list_init <-  lapply(1:10, function(hhs) lapply(HHStypePopulation[[hhs]], function(x) lapply(1:3,function(l) array(0,dim = c(length(infectiousPeriods),l)))))
HHSIarr_list_init <-  lapply(1:10, function(hhs) lapply(HHStypePopulation[[hhs]], function(x) lapply(1:4,function(i) array(0,dim = c(length(latentPeriods),i)))))

Ntot <- sum(populationData[,2:5])
NHHS <- sapply(1:10,function(hhs) sum(populationData[hhs,]))

for (sim in 1:nsim){

  HHSseedInf <- lapply(1:10, function(x) rmultinom(1,seedInf[x],evec)) ## number infected per age group  at beginning

  HHS_ag_seedinit <- list()
  for (hhs in 1:10){
    pmat <- HHSpmat[[hhs]]
    ag_seedinit <- list()
    for (ag in 1:4){
      ag_seedinit[[ag]] <- array(rmultinom(1,HHSseedInf[[hhs]][ag],pmat[[ag]]),dim = c(length(infectiousPeriods),length(latentPeriods))) 
    }
    HHS_ag_seedinit[[hhs]] <- ag_seedinit
  }
  
  for (hhs in 1:10){
    Earr_list_init <- HHSEarr_list_init[[hhs]]
    ag_seedinit <- HHS_ag_seedinit[[hhs]]
    for (ag in 1:4){
      for (l in 1:3){
        Earr_list_init[[ag]][[l]][,1] <- ag_seedinit[[ag]][,l]
      }
    }
    HHSEarr_list_init[[hhs]] <- Earr_list_init
  }
  
  parameterList[[sim]] <- list(R0)
  
  HHSIarr <- NULL
  HHSSarr <- NULL
  HHSincarr <- NULL
  
  for (hhs in 1:10){
    
    typePopulation <- HHStypePopulation[[hhs]]
    
    Earr_list_init <- HHSEarr_list_init[[hhs]]
    Iarr_list_init <- HHSIarr_list_init[[hhs]]
    
    E_init <- list()
    for (ag in 1:4){
      mat <- NULL
      for (l in 1:3){
        mat <- cbind(mat,rowSums(Earr_list_init[[ag]][[l]]))
      }
      E_init[[ag]] <- mat
    }
    
    I_init <- list()
    for (ag in 1:4){
      mat <- NULL
      for (i in 1:4){
        mat <- rbind(mat,rowSums(Iarr_list_init[[ag]][[i]]))
      }
      I_init[[ag]] <- mat
    }
    
    Sarr_init <- lapply(1:4,function(ag) typePopulation[[ag]] - E_init[[ag]] - I_init[[ag]])
    HHSSarr[[hhs]] <- Sarr_init
    HHSEarr_list[[hhs]] <- HHSEarr_list_init[[hhs]]
    HHSIarr_list[[hhs]] <- HHSIarr_list_init[[hhs]]
  }
  
  time <- 1
  Itotls <- incls <- NULL
  Sproptotls <- NULL
  newinf_list <- list()
  while (time <= durEpidemic){
    
    for (hhs in 1:10){
      
      outlist <- newInfect(t=time,hhs,HHSIarr_list,HHSSarr,HHSEarr_list,HHStypePopulation,R0,infectiousPeriodMean,HHSfact)
      HHSSarr[[hhs]] <- outlist[[1]]
      HHSEarr_list[[hhs]] <- outlist[[2]]
      newinf_list[[hhs]] <- outlist[[3]]
      
      HHSIarr_list[[hhs]] <- dayProgI(hhs,HHSEarr_list,HHSIarr_list)
      HHSEarr_list[[hhs]] <- dayProgE(hhs,HHSEarr_list)
    }
    time <- time + 1
    incls <- c(incls,sum(unlist(newinf_list)))
    Itotls <- c(Itotls,sum(unlist(HHSIarr_list)))
  }
  inctotarr <- rbind(inctotarr,incls,deparse.level = 0)
  Itotarr <- rbind(Itotarr,Itotls, deparse.level = 0)
}


sink('ParameterList3.txt')
print(parameterList)
sink()

write.csv(t(inctotarr),'Daily incidence in 10 sims, const parms, dependent HHS.csv')

maxy <- max(inctotarr)

xlim <- 140

colls <- rainbow(nsim)

plot(inctotarr[1,1:xlim],type = 'l',col = colls[1],lwd = 2,ylim = c(0,maxy),xlab = 'Outbreak Day',ylab = 'Number Infectious')

for (sim in 2:nsim){
  lines(inctotarr[sim,1:xlim],col = colls[sim], lwd = 2) 
}
