###########################
# SMALLPOX MODEL
###########################

# SEIR <- function(t, x, parms)
# Use: calculates the derivatives for the SEIR model
# Input: 
#      t: time (not used here, because there is no explicit time dependence)
#      x: vector of the current values of all variables 
#      parms: vector of all parameter values
# Output:
#      der: vector of derivatives

# To use the lsoda function, the model function has to be a function of t (time),
# x (the vector of the variables) and parms (the parameters of the model).

SEIR <- function(t, x, parms){
  
  with(as.list(c(parms,x)),{
    Inv <-   ifelse(t<22,Inv,210000)
    Inv <-   ifelse(t<29,Inv,31189)
    Inv <-   ifelse(t<30,Inv,0)
    Inv <-   ifelse(t<33,Inv,210000)
    Inv <-   ifelse(t<43,Inv,151784)
    Inv <-   ifelse(t<44,Inv,0)
    dS2 <- - b*S2*I3
    dS3 <- - b*S3*I3-Inv*S3/(S3+I1+I2+I3)
    dI1 <-   b*S3*I3+b*S2*I3-r1*I1-Inv*I1/(S3+I1+I2+I3)
    dI2 <-   r1*I1-r2*I2
    dI3 <-   r2*I2-r3*I3
    dI4 <-   r3*I3-r4*I4
    dR  <-   r4*I4+Inv*S3/(S3+I1+I2+I3)+Inv*I1/(S3+I1+I2+I3)
    der <-   c(dS2, dS3, dI1, dI2, dI3, dI4,dR)
    list(der)  # the output must be returned           
  }) # end of 'with'
}  # end of function definition

###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
library(deSolve)

### INITIALIZE PARAMETER SETTINGS

parms <- c(b=1/3752973, r1=1/3, r2=1/8, r3=1/3, r4=1/12, Inv=0)    # set the parameters of the model
dt    <- seq(0,120,0.1)      # set the time points for evaluation


# Initial conditions

inits <- c(S2=34818, S3=3717655, I1=500, I2=0, I3=0, I4=0,R=0)     

# Calculation of the total number of individuals in each sub-population

N <- sum(inits)
N
### SIMULATION OF THE MODEL

## Use lsoda to solve the differential equations numerically. 

out <- as.data.frame(lsoda(inits, dt, SEIR, parms=parms)) # this way our set 'parms' will be used as default

Infected_22<-c(out$I2)
Effective_Rep_Number_22<-c((3/3752973)*(cbind(out$S2+out$S3)))



###########################
# SMALLPOX MODEL
###########################

# SEIR <- function(t, x, parms)
# Use: calculates the derivatives for the SEIR model
# Input: 
#      t: time (not used here, because there is no explicit time dependence)
#      x: vector of the current values of all variables 
#      parms: vector of all parameter values
# Output:
#      der: vector of derivatives

# To use the lsoda function, the model function has to be a function of t (time),
# x (the vector of the variables) and parms (the parameters of the model).

SEIR <- function(t, x, parms){
  
  with(as.list(c(parms,x)),{
    Inv <-   ifelse(t<29,Inv,210000)
    Inv <-   ifelse(t<36,Inv,31189)
    Inv <-   ifelse(t<37,Inv,0)
    Inv <-   ifelse(t<40,Inv,210000)
    Inv <-   ifelse(t<50,Inv,151784)
    Inv <-   ifelse(t<51,Inv,0)
    dS2 <- - b*S2*I3
    dS3 <- - b*S3*I3-Inv*S3/(S3+I1+I2+I3)
    dI1 <-   b*S3*I3+b*S2*I3-r1*I1-Inv*I1/(S3+I1+I2+I3)
    dI2 <-   r1*I1-r2*I2
    dI3 <-   r2*I2-r3*I3
    dI4 <-   r3*I3-r4*I4
    dR  <-   r4*I4+Inv*S3/(S3+I1+I2+I3)+Inv*I1/(S3+I1+I2+I3)
    der <-   c(dS2, dS3, dI1, dI2, dI3, dI4,dR)
    list(der)  # the output must be returned           
  }) # end of 'with'
}  # end of function definition

###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
library(deSolve)

### INITIALIZE PARAMETER SETTINGS

parms <- c(b=1/3752973, r1=1/3, r2=1/8, r3=1/3, r4=1/12, Inv=0)    # set the parameters of the model
dt    <- seq(0,120,0.1)      # set the time points for evaluation


# Initial conditions

inits <- c(S2=34818, S3=3717655, I1=500, I2=0, I3=0, I4=0,R=0)     

# Calculation of the total number of individuals in each sub-population

N <- sum(inits)
N
### SIMULATION OF THE MODEL

## Use lsoda to solve the differential equations numerically. 

out <- as.data.frame(lsoda(inits, dt, SEIR, parms=parms)) # this way our set 'parms' will be used as default

Infected_29<-c(out$I2)
Effective_Rep_Number_29<-c((3/3752973)*(cbind(out$S2+out$S3)))




###########################
# SMALLPOX MODEL
###########################

# SEIR <- function(t, x, parms)
# Use: calculates the derivatives for the SEIR model
# Input: 
#      t: time (not used here, because there is no explicit time dependence)
#      x: vector of the current values of all variables 
#      parms: vector of all parameter values
# Output:
#      der: vector of derivatives

# To use the lsoda function, the model function has to be a function of t (time),
# x (the vector of the variables) and parms (the parameters of the model).

SEIR <- function(t, x, parms){
  
  with(as.list(c(parms,x)),{
    Inv <-   ifelse(t<36,Inv,210000)
    Inv <-   ifelse(t<43,Inv,31189)
    Inv <-   ifelse(t<44,Inv,0)
    Inv <-   ifelse(t<47,Inv,210000)
    Inv <-   ifelse(t<57,Inv,151784)
    Inv <-   ifelse(t<58,Inv,0)
    dS2 <- - b*S2*I3
    dS3 <- - b*S3*I3-Inv*S3/(S3+I1+I2+I3)
    dI1 <-   b*S3*I3+b*S2*I3-r1*I1-Inv*I1/(S3+I1+I2+I3)
    dI2 <-   r1*I1-r2*I2
    dI3 <-   r2*I2-r3*I3
    dI4 <-   r3*I3-r4*I4
    dR  <-   r4*I4+Inv*S3/(S3+I1+I2+I3)+Inv*I1/(S3+I1+I2+I3)
    der <-   c(dS2, dS3, dI1, dI2, dI3, dI4,dR)
    list(der)  # the output must be returned           
  }) # end of 'with'
}  # end of function definition

###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
library(deSolve)

### INITIALIZE PARAMETER SETTINGS

parms <- c(b=1/3752973, r1=1/3, r2=1/8, r3=1/3, r4=1/12, Inv=0)    # set the parameters of the model
dt    <- seq(0,120,0.1)      # set the time points for evaluation


# Initial conditions

inits <- c(S2=34818, S3=3717655, I1=500, I2=0, I3=0, I4=0,R=0)     

# Calculation of the total number of individuals in each sub-population

N <- sum(inits)
N
### SIMULATION OF THE MODEL

## Use lsoda to solve the differential equations numerically. 

out <- as.data.frame(lsoda(inits, dt, SEIR, parms=parms)) # this way our set 'parms' will be used as default
Infected_36<-c(out$I2)
Effective_Rep_Number_36<-c((3/3752973)*(cbind(out$S2+out$S3)))


###########################
# SMALLPOX MODEL
###########################

# SEIR <- function(t, x, parms)
# Use: calculates the derivatives for the SEIR model
# Input: 
#      t: time (not used here, because there is no explicit time dependence)
#      x: vector of the current values of all variables 
#      parms: vector of all parameter values
# Output:
#      der: vector of derivatives

# To use the lsoda function, the model function has to be a function of t (time),
# x (the vector of the variables) and parms (the parameters of the model).

SEIR <- function(t, x, parms){
  
  with(as.list(c(parms,x)),{
    Inv <-   ifelse(t<43,Inv,210000)
    Inv <-   ifelse(t<50,Inv,31189)
    Inv <-   ifelse(t<51,Inv,0)
    Inv <-   ifelse(t<54,Inv,210000)
    Inv <-   ifelse(t<64,Inv,151784)
    Inv <-   ifelse(t<65,Inv,0)
    dS2 <- - b*S2*I3
    dS3 <- - b*S3*I3-Inv*S3/(S3+I1+I2+I3)
    dI1 <-   b*S3*I3+b*S2*I3-r1*I1-Inv*I1/(S3+I1+I2+I3)
    dI2 <-   r1*I1-r2*I2
    dI3 <-   r2*I2-r3*I3
    dI4 <-   r3*I3-r4*I4
    dR  <-   r4*I4+Inv*S3/(S3+I1+I2+I3)+Inv*I1/(S3+I1+I2+I3)
    der <-   c(dS2, dS3, dI1, dI2, dI3, dI4,dR)
    list(der)  # the output must be returned           
  }) # end of 'with'
}  # end of function definition

###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
library(deSolve)

### INITIALIZE PARAMETER SETTINGS

parms <- c(b=1/3752973, r1=1/3, r2=1/8, r3=1/3, r4=1/12, Inv=0)    # set the parameters of the model
dt    <- seq(0,120,0.1)      # set the time points for evaluation


# Initial conditions

inits <- c(S2=34818, S3=3717655, I1=500, I2=0, I3=0, I4=0,R=0)     

# Calculation of the total number of individuals in each sub-population

N <- sum(inits)
N
### SIMULATION OF THE MODEL

## Use lsoda to solve the differential equations numerically. 

out <- as.data.frame(lsoda(inits, dt, SEIR, parms=parms)) # this way our set 'parms' will be used as default
Infected_43<-c(out$I2)
Effective_Rep_Number_43<-c((3/3752973)*(cbind(out$S2+out$S3)))


###########################
# SMALLPOX MODEL
###########################

# SEIR <- function(t, x, parms)
# Use: calculates the derivatives for the SEIR model
# Input: 
#      t: time (not used here, because there is no explicit time dependence)
#      x: vector of the current values of all variables 
#      parms: vector of all parameter values
# Output:
#      der: vector of derivatives

# To use the lsoda function, the model function has to be a function of t (time),
# x (the vector of the variables) and parms (the parameters of the model).

SEIR <- function(t, x, parms){
  
  with(as.list(c(parms,x)),{
    Inv <-   ifelse(t<50,Inv,277300)
    Inv <-   ifelse(t<57,Inv,31189)
    Inv <-   ifelse(t<58,Inv,0)
    Inv <-   ifelse(t<61,Inv,210000)
    Inv <-   ifelse(t<71,Inv,151784)
    Inv <-   ifelse(t<72,Inv,0)
    dS2 <- - b*S2*I3
    dS3 <- - b*S3*I3-Inv*S3/(S3+I1+I2+I3)
    dI1 <-   b*S3*I3+b*S2*I3-r1*I1-Inv*I1/(S3+I1+I2+I3)
    dI2 <-   r1*I1-r2*I2
    dI3 <-   r2*I2-r3*I3
    dI4 <-   r3*I3-r4*I4
    dR  <-   r4*I4+Inv*S3/(S3+I1+I2+I3)+Inv*I1/(S3+I1+I2+I3)
    der <-   c(dS2, dS3, dI1, dI2, dI3, dI4,dR)
    list(der)  # the output must be returned           
  }) # end of 'with'
}  # end of function definition

###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
library(deSolve)

### INITIALIZE PARAMETER SETTINGS

parms <- c(b=1/3752973, r1=1/3, r2=1/8, r3=1/3, r4=1/12, Inv=0)    # set the parameters of the model
dt    <- seq(0,120,0.1)      # set the time points for evaluation


# Initial conditions

inits <- c(S2=34818, S3=3717655, I1=500, I2=0, I3=0, I4=0,R=0)     

# Calculation of the total number of individuals in each sub-population

N <- sum(inits)
N
### SIMULATION OF THE MODEL

## Use lsoda to solve the differential equations numerically. 

out <- as.data.frame(lsoda(inits, dt, SEIR, parms=parms)) # this way our set 'parms' will be used as default
Infected_50<-c(out$I2)
Effective_Rep_Number_50<-c((3/3752973)*(cbind(out$S2+out$S3)))

##Calculating the vector of total infected 
Total_Infected<-cbind(dt,Infected_22,Infected_29,Infected_36,Infected_43,Infected_50)

##Calculating the vector of the various Effective Reproduction Numbers 

Effective_Rep_Number_all<-cbind(dt,Effective_Rep_Number_22,Effective_Rep_Number_29,Effective_Rep_Number_36,Effective_Rep_Number_43,Effective_Rep_Number_50)

library(xlsx)

write.xlsx(Total_Infected, "C:/Users/oresths/Desktop/Total_Infected.xlsx")

##Plotting the number of infected individuals at stage 2 of infection

matplot(x = Total_Infected[,1], y = Total_Infected[,-1], type = "l", lwd = 2,
        lty = "solid", col = c("black","yellow", "blue", "green","red"),
        xlab = "time (days)", ylab = "Infected individuals", main = "Fourth set of scenarios")

legend("topright", col = c("black","yellow", "blue", "green","red"),
       legend = c("Vaccination at day 22", "Vaccination at day 29", "Vaccination at day 36", "Vaccination at day 43", "Vaccination at day 50"), lwd = 2)
