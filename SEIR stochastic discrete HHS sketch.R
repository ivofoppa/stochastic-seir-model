
HHSseedInf <- lapply(seq_along(HHS_order), function(x) rmultinom(1,seedInf[x],evec)) ## number infected per age group  at beginning
names(HHSseedInf) <- sapply(HHS_order,function(x) toString(x))

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
      for (l in seq_along(latentPeriods)){
        Iarr_list2[[ag]][[i]][l,1] <- Earr_list[[ag]][[l]][i,l]
      }
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

HHSseedInf <- lapply(HHS_order, function(x) rmultinom(1,seedInf[x],evec)) ## number infected per age group  at beginning
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
# Has the numbers of individuals by "type" by day since infection
HHSEarr_list_init <-  lapply(1:10, function(hhs) lapply(typePopulation, function(x) lapply(seq_along(latentPeriods),function(x) array(0,dim = c(length(infectiousPeriods),x)))))

HHSIarr_list_init <-  lapply(1:10, function(hhs) lapply(typePopulation, function(x) lapply(seq_along(infectiousPeriods),function(x) array(0,dim = c(length(latentPeriods),x)))))

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

corrls <- round(runif(10,34-30,34+30)) ## Random peak R0 between ~ Dec 15 and Feb 15


Sarr_init <- lapply(1:4,function(x) typePopulation[[x]] - E_init[[x]] - I_init[[x]])

inits <- list()

parameters <- list(HHSEarr_list_init=HHSEarr_list_init,HHSIarr_list_init=HHSIarr_list_init,
                   corrls=corrls,R0minls=R0minls,R0maxls=R0maxls,HHStypePopulation=HHStypePopulation,
                   HHS_ag_seedinit=HHS_ag_seedinit)

hhs <- 2

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

while (time < durEpidemic){
  outlist <- newInfect(t=time,Iarr_list,Sarr,Earr_list,typePopulation,R0min,R0max,corr)
  Sarr <- outlist[[1]]
  Earr_list <- outlist[[2]]
  
  Iarr_list <- dayProgI(Earr_list,Iarr_list)
  Earr_list <- dayProgE(Earr_list)
  Ils <- c(Ils,sum(unlist(Iarr_list)))
  
  time <- time + 1
}
