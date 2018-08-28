# Model parameters

populationFractions <- c(0.06194508, 0.17698821, 0.59678721, 0.16427950) #UK 2011 data, age group ranges: 0-4, 5-19, 20-64, 65+
populationLabels <- c("Ages 0-4", "Ages 5-19", "Ages 20-64", "Ages 65+")
# This is the POLYMOD matrix for UK all contacts, regrouped
PolyMod_matrix <- matrix(nrow = 4, ncol = 4)
PolyMod_matrix[1,] <- c(1.9200000, 0.5203550, 0.4761399, 0.1272683)
PolyMod_matrix[2,] <- c(1.4867476, 8.7525691, 1.7676622, 0.7706686)
PolyMod_matrix[3,] <- c(4.5871960, 5.9603867, 7.8770903, 3.5312914)
PolyMod_matrix[4,] <- c(0.3375179, 0.7153304, 0.9720697, 1.8867659)

# Matrix entry (row i, column j) denotes the number of potentially infectious contacts a single individual
# from group j has with individuals from group i each day. This will be re-scaled to have spectral radius 1.
contactMatrix <- PolyMod_matrix

R0_High <- 1.65 # Gives a serologic AR of ~60%; symptomatic AR of ~30%
R0_Low <- 1.35 # Gives a serologic AR of ~40%; symptomatic AR of ~20%
latentPeriod <- 1.0 # Days
infectiousPeriod <- 1.5 # Days

seedInfections <- 100

baselineModel <- SEIRModel(population = 310e6,
                           populationFractions = populationFractions,
                           contactMatrix = contactMatrix,
                           R0 = R0_High,
                           latentPeriod = latentPeriod,
                           infectiousPeriod = infectiousPeriod,
                           seedInfections = seedInfections)

# Type these in the command line to see interactive information
# baselineModel
# plot(baselineModel, incidence = TRUE, populationLabels = populationLabels)

fractionSymptomatic <- 0.5

caseHospitalizationRatio_Low <- 1.2 / 100 # 'Low severity' 0.8-1.5% (Used midpoint)
caseFatalityRatio_Low <- 0.08 / 100 # 'Low severity' 0.05-0.1% (Used midpoint)

caseHospitalizationRatio_High <- 4 / 100 # 'High severity' 3-5% (Used midpoint)
caseFatalityRatio_High <- 0.38 / 100 # 'High severity' 0.25-0.5% (Used midpoint)

baselineHospitalizationsByAge <- getHospitalizations(baselineModel,
                                                     fractionSymptomatic = fractionSymptomatic,
                                                     caseHospitalizationRatio = caseHospitalizationRatio_High)
baselineHospitalizationsOverall <- getHospitalizations(baselineModel,
                                                       fractionSymptomatic = fractionSymptomatic,
                                                       caseHospitalizationRatio = caseHospitalizationRatio_High,
                                                       byGroup = FALSE)

baselineDeathsByAge <- getDeaths(baselineModel,
                                 fractionSymptomatic = fractionSymptomatic,
                                 caseFatalityRatio = caseFatalityRatio_High)
baselineDeathsOverall <- getDeaths(baselineModel,
                                   fractionSymptomatic = fractionSymptomatic,
                                   caseFatalityRatio = caseFatalityRatio_High,
                                   byGroup = FALSE)

# Vaccine Modeling

vaccineDosesAvailable <- 310e6 * 0.7 * 2 # Vaccine for 70% of population, 2 doses each
vaccineCampaignStartDay <- 2 * 7 # Start after 2 weeks
vaccineAvailabilityByDay <- c(rep(0, vaccineCampaignStartDay), vaccineDosesAvailable)
vaccinationRatePerWeek_Low <- 10e6
vaccinationRatePerWeek_High <- 20e6

vaccineEfficacyDose1 <- 0.0 # First dose does nothing
vaccineEfficacyDose2 <- c(0.62, 0.62, 0.62, 0.43) # 62% efficacy for everyone below age 60, and 43% for those above
vaccineEfficacyDelay <- 7 
dose2Delay <- 14

# By default, vaccine uptake is proportional to population
# Child (Ages 0-19) fraction of population is:
childPopulationFraction <- sum(populationFractions[1:2])
# So, to have 50% of vaccine taken by children, adjust uptake:
vaccineUptakeMultiplierByAge_5050 <- c((1 - childPopulationFraction) / childPopulationFraction,
                                       (1 - childPopulationFraction) / childPopulationFraction,
                                       1,
                                       1)
# To have 10% of vaccine taken by children, adjust again:
vaccineUptakeMultiplierByAge_1090 <- c((1 - childPopulationFraction) / childPopulationFraction,
                                       (1 - childPopulationFraction) / childPopulationFraction,
                                       9,
                                       9)

vaccineModel <- SEIRV2DoseModel(population = 310e6,
                                populationFractions = populationFractions,
                                contactMatrix = contactMatrix,
                                R0 = R0_High,
                                latentPeriod = latentPeriod,
                                infectiousPeriod = infectiousPeriod,
                                seedInfections = seedInfections,
                                vaccineAvailabilityByDay = vaccineAvailabilityByDay,
                                vaccineAdministrationRatePerDay = vaccinationRatePerWeek_High / 7,
                                vaccineUptakeMultiplier = vaccineUptakeMultiplierByAge_5050,
                                dose2Delay = dose2Delay,
                                VEs1 = vaccineEfficacyDose1,
                                VEs2 = vaccineEfficacyDose2,
                                vaccineEfficacyDelay = vaccineEfficacyDelay)

# Results

baselineClinicalCases <- getInfections(baselineModel,
                                       symptomatic = TRUE,
                                       fractionSymptomatic = fractionSymptomatic,
                                       byGroup = FALSE)
clinicalCasesAverted <- baselineClinicalCases - getInfections(vaccineModel,
                                                              symptomatic = TRUE,
                                                              fractionSymptomatic = fractionSymptomatic,
                                                              byGroup = FALSE)
hospitalizationsAverted <- baselineHospitalizationsOverall - getHospitalizations(vaccineModel,
                                                                                 fractionSymptomatic = fractionSymptomatic,
                                                                                 caseHospitalizationRatio = caseHospitalizationRatio_High,
                                                                                 byGroup = FALSE)
deathsAverted <- baselineDeathsOverall - getDeaths(vaccineModel,
                                                   fractionSymptomatic = fractionSymptomatic,
                                                   caseFatalityRatio = caseFatalityRatio_High,
                                                   byGroup = FALSE)

par(mfrow = c(1, 2))
plot(baselineModel, incidence = TRUE, populationLabels = populationLabels)
plot(vaccineModel, incidence = TRUE, populationLabels = populationLabels)
par(mfrow = c(1, 1))