
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
## Construct population:
populationData <- read.table('pop_reg5.dat',header = T)

## Creating population according to latent (rows) and infectious periods (columns) in each of the age groups
r0 <- function(t,R0min,R0max,corr){
  y <- (sin((t - corr)/365 * 2 * pi ) + 1) / 2 * (R0max - R0min)
  return(y + R0min)
}

beta <- function(t,R0min,R0max,corr) {
  r0(t,R0min,R0max,corr) / infectiousPeriod / specrad
}


pSE <- function(t,Iarr_list,typePopulation,R0min,R0max,corr){ ## k is age group, I is matrix of # infectious by latent/infectious time type
  I <- sapply(1:4, function(x) sum(unlist(Iarr_list[[x]])))
  pop <- sapply(1:4, function(x) sum(typePopulation[[x]]))
  lambda <- beta(t,R0min,R0max,corr) * contactMatrix %*% as.numeric(I/pop)
  return (as.numeric(1-exp(-lambda)))
}

newInfect <- function(t,Iarr_list,Sarr,Earr_list,typePopulation,R0min,R0max,corr){
  newinfList <- list()
  Sarr2 <- Sarr
  Earr_list2 <- Earr_list
  pinfls <- pSE(t,Iarr_list,typePopulation,R0min,R0max,corr)
  for (ag in 1:4){
    Sarr_ag <- Sarr[[ag]]
    Earr <- Earr_list2[[ag]]
    pinf <- pinfls[[ag]]
    newinf <- structure(vapply(Sarr_ag, function(x) rbinom(1,x,pinf), numeric(1)), dim=dim(Sarr_ag))
    Earr_list2[[ag]] <- sapply(seq_along(Earr), function(x) {Earr[[x]][,1] <- Earr[[x]][,1] + newinf[,x]; return(Earr[[x]])})
    Sarr2[[ag]] <- Sarr2[[ag]] - newinf
    newinfList[[ag]] <- newinf
  }
  return(list(Sarr2,Earr_list2,newinfList))
}

dayProgI <- function(Earr_list,Iarr_list){
  Iarr_list2 <- Iarr_list
  for (ag in 1:4){
    for (i in seq_along(infectiousPeriods)){
      newinf <- sapply(seq_along(latentPeriods),function(l) Earr_list[[ag]][[l]][i,l])
      Iarr_list2[[ag]][[i]] <- cbind(newinf,Iarr_list2[[ag]][[i]][,1:i-1],deparse.level = 0)
      }
  }
  return(Iarr_list2)
} 

dayProgE <- function(Earr_list){
  Earr_list2 <- Earr_list
  for(ag in 1:4){
    for (l in seq_along(latentPeriods)){
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
seedInf <- rep(10,10)
seedInf[c(2,6)] <- 500
seedInf[c(3,4)] <- 50

nsim <- 10
colls <- rainbow(nsim)

# Has the numbers of individuals by "type" by day since infection
for (sim in 1:nsim){
  HHSEarr_list_init <-  lapply(1:10, function(hhs) lapply(typePopulation, function(x) lapply(seq_along(latentPeriods),function(l) array(0,dim = c(length(infectiousPeriods),l)))))
  HHSIarr_list_init <-  lapply(1:10, function(hhs) lapply(typePopulation, function(x) lapply(seq_along(infectiousPeriods),function(i) array(0,dim = c(length(latentPeriods),i)))))
  
  HHSseedInf <- lapply(1:10, function(x) rmultinom(1,seedInf[x],evec)) ## number infected per age group  at beginning
  names(HHSseedInf) <- sapply(1:10,function(x) toString(x))
  
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
      for (l in seq_along(latentPeriods)){
        larr <- length(Earr_list_init[[ag]][[l]][1,])
        Earr_list_init[[ag]][[l]][,larr] <- ag_seedinit[[ag]][,l]
      }
    }
    HHSEarr_list_init[[hhs]] <- Earr_list_init
  }
  
  for (hhs in 1:10){
    Iarr_list_init <- HHSIarr_list_init[[hhs]]
    ag_seedinit <- HHS_ag_seedinit[[hhs]]
    for (ag in 1:4){
      for (i in seq_along(infectiousPeriods)){
        Iarr_list_init[[ag]][[i]][,1] <- ag_seedinit[[ag]][i,]
      }
    }
    HHSIarr_list_init[[hhs]] <- Iarr_list_init
  }
  
  
  
  
  E_init <- list()
  for (ag in 1:4){
    mat <- NULL
    for (l in seq_along(latentPeriods)){
      mat <- cbind(mat,rowSums(Earr_list_init[[ag]][[l]]))
    }
    E_init[[ag]] <- mat
  }
  
  I_init <- list()
  for (ag in 1:4){
    mat <- NULL
    for (i in seq_along(infectiousPeriods)){
      mat <- rbind(mat,rowSums(Iarr_list_init[[ag]][[i]]))
    }
    I_init[[ag]] <- mat
  }
  
  R0minls <- runif(10,1.5,2.2)
  R0maxls <- sapply(R0minls, function(x) runif(1,x,2.2))
  
  corrls <- round(runif(10,64-30,64+30)) ## Random peak R0 between ~ Dec 15 and Feb 15
  
  
  Sarr_init <- lapply(1:4,function(x) typePopulation[[x]] - E_init[[x]] - I_init[[x]])
  
  parameters <- list(HHSEarr_list_init=HHSEarr_list_init,HHSIarr_list_init=HHSIarr_list_init,
                     corrls=corrls,R0minls=R0minls,R0maxls=R0maxls,HHStypePopulation=HHStypePopulation,
                     HHS_ag_seedinit=HHS_ag_seedinit)
  
  # hhs <- 2
  HHSIarr <- NULL
  for (hhs in 1:10){
    latentPeriod <- runif(1,1.3,1.7)
    
    lp3 <- -1; lp2 <- -1
    while (lp3 < 0 | lp3 > 1 | lp2 < 0 | lp2 > 1 ){
      lp1 <- runif(1,.1,.9)
      lp2 <- (latentPeriod - 3 + 2*lp1)/(-1)
      lp3 <- 1 - lp1 - lp2
    }
    
    latentPeriods <- c(lp1,lp2,lp3) ## frequencies of latent period durations 1-3 days
    
    ## Combinations of infectious/latent periods (discrete) are considered "types"
    infectiousPeriod <- runif(1,2.2,2.8)
    
    ip4 <- -1; ip3 <- -1
    while (ip4 < 0 | ip4 > 1 | ip3 < 0 | ip3 > 1 ){
      ip1 <- runif(1,0,.9)
      ip2 <- runif(1,0,1 - ip1)
      ip3 <- 4 - infectiousPeriod - 3*ip1 - 2*ip2
      ip4 <- 1 - ip1 - ip2 - ip3
    }
    
    infectiousPeriods <- c(ip1,ip2,ip3,ip4) ## frequencies of infectious durations 3-6 days
    
    
    time <- 1
    Ils <- NULL
    
    Earr_list <- parameters$HHSEarr_list_init[[hhs]]
    Iarr_list <- parameters$HHSIarr_list_init[[hhs]]
    
    E_init <- list()
    for (ag in 1:4){
      mat <- NULL
      for (l in seq_along(latentPeriods)){
        mat <- cbind(mat,rowSums(Earr_list[[ag]][[l]]))
      }
      E_init[[ag]] <- mat
    }
    
    I_init <- list()
    for (ag in 1:4){
      mat <- NULL
      for (i in seq_along(infectiousPeriods)){
        mat <- rbind(mat,rowSums(Iarr_list[[ag]][[i]]))
      }
      I_init[[ag]] <- mat
    }
    
    typePopulation <- parameters$HHStypePopulation[[hhs]]
    Sarr <- lapply(1:4, function(x) typePopulation[[x]] - I_init[[x]] - E_init[[x]]) 
    
    R0min <- R0minls[hhs]
    R0max <- R0maxls[hhs]
    corr <- corrls[hhs]
    
    while (time <= durEpidemic){
      outlist <- newInfect(t=time,Iarr_list,Sarr,Earr_list,typePopulation,R0min,R0max,corr)
      Sarr <- outlist[[1]]
      Earr_list <- outlist[[2]]
      
      Iarr_list <- dayProgI(Earr_list,Iarr_list)
      Earr_list <- dayProgE(Earr_list)
      Ils <- c(Ils,sum(unlist(Iarr_list)))
      
      time <- time + 1
    }
    HHSIarr <- rbind(HHSIarr,Ils,deparse.level = 0)
  }
  
  if (sim==1) {plot(colSums(HHSIarr),type = 'l',col = colls[sim])}
  else lines(colSums(HHSIarr),col = colls[sim])
}
