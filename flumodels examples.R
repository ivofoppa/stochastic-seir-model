library(flumodels)

model <- SEIRModel(R0 = 2.2,
                   latentPeriod = 1.5,
                   infectiousPeriod = 2.5,
                   seedInfections = 0.0001 * 310e6,
                   simulationLength = 300,
                   population = 310e6)

model2 <- SEIRModel(R0 = 1.3,
                    latentPeriod = 1.5,
                    infectiousPeriod = 2.5,
                    seedInfections = 0.0001 * 310e6,
                    simulationLength = 300,
                    population = 310e6,
                    populationFractions = c(0.25, 0.75))

model3 <- SEIRModel(R0 = 1.3,
                    latentPeriod = 1.5,
                    infectiousPeriod = 2.5,
                    seedInfections = 0.0001 * 310e6,
                    simulationLength = 300,
                    population = 310e6,
                    populationFractions = c(0.25, 0.75),
                    contactMatrix = matrix(c(18,  3,
                                              9, 12), ncol=2, byrow = TRUE))

model4 <- SEIRVModel(R0 = 1.3,
                     latentPeriod = 1.5,
                     infectiousPeriod = 2.5,
                     seedInfections = 0.0001 * 310e6,
                     simulationLength = 300,
                     population = 310e6,
                     populationFractions = c(0.25, 0.75),
                     contactMatrix = matrix(c(18,  3,
                                              9, 12), ncol=2, byrow = TRUE),
                     vaccineAvailabilityByDay = 40e6,
                     vaccineAdministrationRatePerDay = 2e6,
                     VEs = 0.5,
                     VEi = 0.5,
                     vaccineEfficacyDelay = 0)

model5 <- SEIRV2DoseModel(R0 = 1.3,
                          latentPeriod = 1.5,
                          infectiousPeriod = 2.5,
                          seedInfections = 0.0001 * 310e6,
                          simulationLength = 300,
                          population = 310e6,
                          populationFractions = c(0.25, 0.75),
                          contactMatrix = matrix(c(18,  3,
                                                   9, 12), ncol=2, byrow = TRUE),
                          vaccineAvailabilityByDay = 40e6,
                          vaccineAdministrationRatePerDay = 2e6,
                          dose2Delay = 21,
                          VEs2 = 0.5,
                          VEi2 = 0.5,
                          vaccineEfficacyDelay = 0)

model6 <- SEIRTModel(R0 = 1.3,
                     latentPeriod = 1.5,
                     infectiousPeriod = 2.5,
                     seedInfections = 0.0001 * 310e6,
                     simulationLength = 300,
                     population = 310e6,
                     populationFractions = c(0.25, 0.75),
                     contactMatrix = matrix(c(18,  3,
                                              9, 12), ncol=2, byrow = TRUE),
                     fractionSymptomatic = 0.5,
                     fractionSeekCare = 0.3,
                     fractionDiagnosedAndPrescribedOutpatient = 0.5,
                     fractionAdhere = 0.2,
                     fractionAdmitted = 1.0,
                     fractionDiagnosedAndPrescribedInpatient = 0.5,
                     AVEi = 0.5,
                     AVEp = 0.5)

# Plots two models against each other
comparePlot <- function(model1, model2) {
  model1Times <- model1$rawOutput[, "time"]
  model1Infections <- getInfectionTimeSeries(model1, byGroup = FALSE)
  model2Times <- model1$rawOutput[, "time"]
  model2Infections <- getInfectionTimeSeries(model2, byGroup = FALSE)
  
  scalePoints <- pretty(c(0, max(model1Infections, model2Infections, na.rm = TRUE)))
  ylim <- c(0, max(scalePoints))

  plot(model1Times, model1Infections, type = 'l', xlab = "Time (days)", ylab = "",
       ylim = ylim, yaxt = "n", lwd = 3)
  lines(model2Times, model2Infections, lwd = 3, col = 2)
  title(main = "Infection Prevalence")
  
  # Make a legend
  legend("topright", legend = c(deparse(substitute(model1)), deparse(substitute(model2))),
         col = c(1, 2), bty = "n", lwd = 2)
  
  # Scale y-axis
  exponent <- floor(log10(ylim[2]) / 3)
  if (exponent < 2) { #Hundreds of thousands or fewer
    axis(2, at = scalePoints, las = 1, lab = format(scalePoints, big.mark=",", scientific=FALSE))
  } else if (exponent < 3) { #Millions
    axis(2, at = scalePoints, las = 1, lab = scalePoints / 1e6)
    mtext("Millions", line = 3, side = 2, adj = 1)
  } else if (exponent < 4) { #Billions
    axis(2, at = scalePoints, las = 1, lab = scalePoints / 1e9)
    mtext("Billions", line = 3, side = 2, adj = 1)
  } else {
    axis(2, at = scalePoints, las = 1, lab = format(scalePoints, scientific=TRUE))
  }
}


# plot(model)
# plot(model2)
# plot(model3)
# etc
comparePlot(model3, model4)
