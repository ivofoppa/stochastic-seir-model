library(flumodels)

model <- SEIRModel(R0 = 1.5,
                   latentPeriod = 1.5,
                   infectiousPeriod = 2.5,
                   seedInfections = 1000,
                   simulationLength = 300,
                   population = 310e6,populationFractions = populationFractions,contactMatrix = contactMatrix)


modelTimes <- model$rawOutput[, "time"]
modelInfections <- getInfectionTimeSeries(model, byGroup = FALSE)

selind <- which(modelTimes <= 150)
plot(modelTimes[selind],modelInfections[selind],type = 'l')
lines(modelTimes[selind],modelInfections[selind],type = 'l')
