library(flumodels)

model <- SEIRModel(R0 = 2.3,
                   latentPeriod = 1.5,
                   infectiousPeriod = 2.5,
                   seedInfections = 1000,
                   simulationLength = 300,priorImmunity = 0,
                   population = 310e6,populationFractions = populationFractions,contactMatrix = contactMatrix)


modelTimes <- model$rawOutput[, "time"]
modelInfections <- getInfectionTimeSeries(model, byGroup = FALSE)
modelSusceptibles <- getSus
selind <- which(modelTimes <= 150)
plot(modelTimes[selind],modelInfections[selind],type = 'l')
lines(modelTimes[selind],modelInfections[selind],type = 'l',lwd = 2, col = 'blue')
