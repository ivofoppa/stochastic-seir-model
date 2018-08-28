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

latentPeriods <- list(d1=.3,d2=.5,d3=.2) ## frequencies of latent period durations 1-3 days

## Combinations of infectious/latent periods (discrete) are considered "types"
infectiousPeriods <- list(d3=.3,d4=.4,d5=.2,d6=.1) ## frequencies of infectious durations 3-6 days
MeaninfectiousPeriod <- sum(unlist(infectiousPeriods)*c(3:6))   ## Mean infectious period

## Construct population:
Npop <- 310e6
Npopag <- rmultinom(1,Npop,c(0.06194508, 0.17698821, 0.59678721, 0.16427950))

## Creating population according to latent (rows) and infectious periods (columns) in each of the age groups
pop <- sapply(Npopag, function(nag) list(sapply(rmultinom(1,nag,c(.3,.5,.2)),function(x) sapply(x, function(y) rmultinom(1,y,c(.3,.4,.2,.1))))))

pmat <- sapply(seq_along(pop),function(x) list(pop[[x]]/sum(pop[[x]])))
R0 <- 1.5
beta <- R0 / MeaninfectiousPeriod / max(eigen(contactMatrix, 
                                              symmetric = FALSE, only.values = TRUE)$values)

parameters <- list(beta=beta,latentPeriods=latentPeriods,infectiousPeriods=infectiousPeriods,
                   pop=pop,pmat=pmat)

seedInf <- 100 ## number infected at beginning

ag_seed <- 3 ## age group in which infection is seeded

E_ag_seedinit <- array(rmultinom(1,seedInf,pmat[[3]]),dim = c(length(infectiousPeriods),length(latentPeriods))) ## Infection only seeded in the "seed" age group

Einit <- lapply(seq_along(pop), function(x) if(x==ag_seed) E_ag_seedinit else pop[[x]]*0)
Sinit <- lapply(seq_along(pop),function(x) pop[[x]] - Einit[[x]])
Iinit <- lapply(seq_along(pop),function(x) pop[[x]]*0)
Rinit <- lapply(seq_along(pop),function(x) pop[[x]]*0)

inits <- list(Sinit=Sinit,Einit=Einit,Iinit=Iinit,Rinit=Rinit)

Ivec <- function(Iarr){
  unlist(lapply(Iinit, function(x) sum(x)))
}

# Has the numbers of individuals by "type" by day since infection
Earr_list <-  lapply(pop, function(x) lapply(seq_along(latentPeriods),function(x) array(0,dim = c(length(infectiousPeriods),x))))
for (k in seq_along(latentPeriods)){
  Earr_list[[ag_seed]][[k]][,1] <- E_ag_seedinit[,k]
}

# Has the numbers of individuals by "type" by day since infection
Iarr_list <-  lapply(pop, function(x) lapply(seq_along(infectiousPeriods),function(x) array(0,dim = c(length(infectiousPeriods), x + (as.numeric(substring(names(infectiousPeriods)[1],2,2)) - 1)))))

## Function to collect Earr_list into 
pSE <- function(k,Iarr_list){ ## k is age group, I is matrix of # infectious by latent/infectious time type
  I <- sapply(seq_along(Iarr_list), function(x) sum(unlist(Iarr_list[[x]])))
  lambda <- contactMatrix[k,] %*% I/sum(pop[[k]])
  return(as.numeric(1-exp(-lambda)))
}

## Function for creating new infections per "type" - latent stage; in age group k
newInfect <- function(k,Iarr_list,Sarr){
  Sarr_ag <- Sarr[[k]]
  pinf <- pSE(k,Iarr_list)
  return(structure(vapply(Sarr_ag, function(x) rbinom(1,x,pinf), numeric(1)), dim=dim(m)))
}


## Function for the stage progression of the latent stage
flowEI <- function(k,Earr,Iarr){
  Ivec <- sapply(seq_along(I), function(x) sum(I[[x]]))
  lambda <- contactMatrix[k,] %*% Ivec
  return(as.numeric(1-exp(-lambda)))
}

SEIRModel2 <- function(population, populationFractions, contactMatrix, R0,
                      latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                      useCommunityMitigation, communityMitigationStartDay,
                      communityMitigationDuration, communityMitigationMultiplier,
                      simulationLength, seedStartDay, tolerance, method) {
  #Check inputs
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIR", argumentList) #Get parameters from checked inputs
  
  initialState <- with(parameters, {
    c(S  = (1 - priorImmunity) * populationFractions,
      E  = 0 * populationFractions,
      I  = 0 * populationFractions,
      R  = priorImmunity * populationFractions)
  })
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIR,
                              seedFunction = doSeed.SEIR)
  
  #Build the SEIRModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- "SEIRModel"
  return(model)
}

#' @title Check SEIR inputs
#' @description Checks the input parameters for the SEIR model
#' @return List of parameters for the SEIR model
#' @keywords internal
checkInputs.SEIR <- function(population = 1, populationFractions = 1, contactMatrix, R0,
                             latentPeriod, infectiousPeriod, seedInfections = 0, priorImmunity = 0,
                             useCommunityMitigation = FALSE, communityMitigationStartDay, 
                             communityMitigationDuration, communityMitigationMultiplier, 
                             simulationLength = 240, seedStartDay = 0, tolerance = 1e-8, method = "lsoda", ...) {
  #population
  checkPositiveNumber(population)
  #populationFractions
  if (!((abs(sum(populationFractions) - 1) < tolerance) && all(populationFractions >= 0))) {
    stop("populationFractions must be positive and sum to 1.", call. = FALSE)
  }
  populationFractions <- as.vector(populationFractions) #Drop names to prevent errors in internal data storage
  #contactMatrix
  if (missing(contactMatrix)) { #If contact matrix is not supplied, default to proportional mixing
    contactMatrix <- as.matrix(populationFractions)[ , rep(1, length(populationFractions)), drop = FALSE]
  } else {
    if (!(all(contactMatrix >=0) && (ncol(contactMatrix) == nrow(contactMatrix))))  {
      stop("contactMatrix must be square and non-negative.", call. = FALSE)
    }
    if (length(populationFractions) != ncol(contactMatrix)) {
      stop("contactMatrix dimensions do not match populationFractions.", call. = FALSE)
    }
  } 
  #R0
  if (missing(R0)) {
    stop("R0 must be specified.", call. = FALSE)
  }
  checkPositive(R0)
  #latentPeriod
  if (missing(latentPeriod)) {
    stop("latentPeriod must be specified.", call. = FALSE)
  }
  checkPositive(latentPeriod)
  #infectiousPeriod
  if (missing(infectiousPeriod)) {
    stop("infectiousPeriod must be specified.", call. = FALSE)
  }
  checkPositive(infectiousPeriod)
  #seedInfections
  checkNonNegative(seedInfections)
  checkDimensionsMatch(seedInfections, populationFractions)
  if (length(seedInfections) == 1) {
    if (seedInfections > population) {
      stop("seedInfections can not exceed the population.", call. = FALSE)
    }
    seedInfections <- seedInfections * populationFractions #Distribute seed infections among population groups proportionately
  } else if (!all(seedInfections <= (population * populationFractions))) {
    stop("seedInfections can not exceed the population by group.", call. = FALSE)
  }
  #priorImmunity
  checkBetween0and1(priorImmunity)
  checkDimensionsMatch(priorImmunity, populationFractions)
  #Community Mitigation
  if (useCommunityMitigation) {
    if (missing(communityMitigationStartDay)) {
      stop("communityMitigationStartDay must be specified when using community mitigation.", 
           call. = FALSE)
    }
    checkNonNegative(communityMitigationStartDay)
    if (missing(communityMitigationDuration)) {
      stop("communityMitigationDuration must be specified when using community mitigation.", 
           call. = FALSE)
    }
    checkNonNegative(communityMitigationStartDay)
    if (missing(communityMitigationMultiplier)) {
      stop("communityMitigationMultiplier must be specified when using community mitigation.", 
           call. = FALSE)
    }
    if (!all(dim(communityMitigationMultiplier) == dim(contactMatrix))) {
      stop("Dimensions of communityMitigationMultiplier do not match those of contactMatrix", 
           call. = FALSE)
    }
    checkNonNegative(communityMitigationMultiplier)
  }
  
  #Collect and return the parameters
  parameters <- list(population = population,
                     populationFractions = populationFractions,
                     contactMatrix = contactMatrix,
                     beta = R0 / infectiousPeriod / max(Mod(eigen(contactMatrix, symmetric = FALSE, only.values = TRUE)$values)),
                     gamma = 1 / infectiousPeriod,
                     lambda = 1 / latentPeriod,
                     seedInfections = seedInfections,
                     priorImmunity = priorImmunity,
                     useCommunityMitigation = useCommunityMitigation,
                     simulationLength = simulationLength,
                     seedStartDay = seedStartDay,
                     tolerance = tolerance,
                     method = method)
  if (useCommunityMitigation) {
    parameters <- append(parameters,
                         list(communityMitigationStartDay = communityMitigationStartDay, 
                              communityMitigationEndDay = communityMitigationStartDay + communityMitigationDuration,
                              communityMitigationMultiplier = communityMitigationMultiplier))
  }
  return(parameters)
}

#This is a utility function that reconstructs the model state as a list so that equations can refer to compartments by name
reconstructState.SEIR <- function(state) {
  numberOfClasses <- length(state)/4 #Each of the 4 classes are vectors of the same length
  S  <- state[                       1 :     numberOfClasses ]
  E  <- state[    (numberOfClasses + 1):(2 * numberOfClasses)]
  I  <- state[(2 * numberOfClasses + 1):(3 * numberOfClasses)]
  R  <- state[(3 * numberOfClasses + 1):(4 * numberOfClasses)]
  return(as.list(environment()))
}

#This function implements the multivariate derivative of the SEIR model
#parameters should define populationFractions, normalizedContactMatrix, beta, lambda, and gamma
#Note that the total population is normalized to be 1
getDerivative.SEIR <- function(t, state, parameters) {
  stateList <- reconstructState.SEIR(state)
  with(append(stateList, parameters), {
    if (useCommunityMitigation) {
      if ((t >= communityMitigationStartDay) && (t < communityMitigationEndDay)) {
        contactMatrix <- communityMitigationMultiplier * contactMatrix
      } 
    }
    #Flows
    forceOfInfection <- beta / populationFractions * (contactMatrix %*% I)
    
    S_to_E <- S * forceOfInfection
    E_to_I <- lambda * E
    I_to_R <- gamma * I
    
    #Derivatives
    dS <- -S_to_E
    dE <- S_to_E - E_to_I
    dI <- E_to_I - I_to_R
    dR <- I_to_R
    
    #Return derivative
    return(list(c(dS, dE, dI, dR)))
  })
}

#This function implements seeding infections in the SEIR model
#parameters should define seedInfections, lambda, and gamma
#Note that the total population is normalized to be 1
doSeed.SEIR <- function(state, parameters) {
  stateList <- reconstructState.SEIR(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions
    E <- E + seedInfectionsFractions / (1 + lambda / gamma)
    I <- I + seedInfectionsFractions / (1 + gamma / lambda)
    #Return derivative
    return(c(S, E, I, R))
  })
}
? 2018 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
Press h to open a hovercard with more details.