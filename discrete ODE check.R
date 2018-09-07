beta <- betafn(1,R0,infectiousPeriod)
parameters <- c(gamma=gamma,delta=delta,beta=beta)
### Define model
### Initial conditions
S1 <- (as.integer(inits$Sinit))[1]
E1 <- (as.integer(inits$Einit))[1]
I1 <- (as.integer(inits$Iinit))[1]
R1 <- (as.integer(inits$Rinit))[1]

S2 <- (as.integer(inits$Sinit))[2]
E2 <- (as.integer(inits$Einit))[2]
I2 <- (as.integer(inits$Iinit))[2]
R2 <- (as.integer(inits$Rinit))[2]

S3 <- (as.integer(inits$Sinit))[3]
E3 <- (as.integer(inits$Einit))[3]
I3 <- (as.integer(inits$Iinit))[3]
R3 <- (as.integer(inits$Rinit))[3]

S4 <- (as.integer(inits$Sinit))[4]
E4 <- (as.integer(inits$Einit))[4]
I4 <- (as.integer(inits$Iinit))[4]
R4 <- (as.integer(inits$Rinit))[4]
time <- 0
# deltat <- .1
I2ls <- NULL
while (time <= durEpidemic){
    # rate of change
    forceOfInfection <- beta / pop * (contactMatrix %*% c(I1,I2,I3,I4))
    
    S_to_E1 <- S1 * forceOfInfection[1]*deltat
    E_to_I1 <- gamma * E1*deltat
    I_to_R1 <- delta * I1*deltat
    
    S_to_E2 <- S2 * forceOfInfection[2]*deltat
    E_to_I2 <- gamma * E2*deltat
    I_to_R2 <- delta * I2*deltat
    
    S_to_E3 <- S3 * forceOfInfection[3]*deltat
    E_to_I3 <- gamma * E3*deltat
    I_to_R3 <- delta * I3*deltat
    
    S_to_E4 <- S4 * forceOfInfection[4]*deltat
    E_to_I4 <- gamma * E4*deltat
    I_to_R4 <- delta * I4*deltat
    
    #erivatives
    S1 <- S1 -S_to_E1
    E1 <- E1 + S_to_E1 - E_to_I1
    I1 <- I1 + E_to_I1 - I_to_R1
    R1 <- R1 + I_to_R1
    
    S2 <- S2 -S_to_E2
    E2 <- E2 + S_to_E2 - E_to_I2
    I2 <- I2 + E_to_I2 - I_to_R2
    R2 <- R2 + I_to_R2
    
    S3 <- S3 -S_to_E3
    E3 <- E3 + S_to_E3 - E_to_I3
    I3 <- I3 + E_to_I3 - I_to_R3
    R3 <- R3 + I_to_R3
    
    S4 <- S4 -S_to_E4
    E4 <- E4 + S_to_E4 - E_to_I4
    I4 <- I4 + E_to_I4 - I_to_R4
    R4 <- R4 + I_to_R4
    
    I2ls <- c(I2ls, I1 + I2 + I3 + I4)
    
    time <- time + deltat
    # return the rate of change
}

toTime <- 150
selind <- which(times <= toTime)

ymax <- max(c(Iarr,I2ls))

colls <- rainbow(nsim)

plot(seq(0,durEpidemic,deltat),I2ls,xlab = 'Outbreak Day',ylab = 'Number Infectious',type = 'l', col = 'blue',xlim = c(0,toTime),ylim = c(0,ymax))
for (sim in 1:10){
  lines(seq(deltat,toTime,deltat),Iarr[sim,seq(deltat,toTime,deltat)/deltat],lwd = 1, col = colls[sim])
}


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
