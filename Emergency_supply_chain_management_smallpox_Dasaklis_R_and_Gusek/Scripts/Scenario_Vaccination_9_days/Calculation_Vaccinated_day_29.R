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
    Inv <-   ifelse(t<29,Inv,277071)
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
dt    <- seq(0,41,0.001)      # set the time points for evaluation


# Initial conditions

inits <- c(S2=34818, S3=3717655, I1=500, I2=0, I3=0, I4=0,R=0)     

# Calculation of the total number of individuals in each sub-population

N <- sum(inits)
N
### SIMULATION OF THE MODEL

## Use lsoda to solve the differential equations numerically. 

out <- as.data.frame(lsoda(inits, dt, SEIR, parms=parms)) # this way our set 'parms' will be used as default
#q<-cbind(out$S2+out$S3)
#r<-c((3/3827124)*(cbind(out$S2+out$S3)))
Effective_Rep_Number<-c((3/3752973)*(cbind(out$S2+out$S3)))

plot(Effective_Rep_Number)

o<-cbind(dt, Effective_Rep_Number)
# Calculate and print the Effective Reproduction Number on the screen

# WRITE THE OUTPUT OF THE SIMULATION IN EXTERNAL ARCHIVES

library(xlsx)
write.xlsx(o, "C:/Users/oresths/Desktop/Effective_Rep_Number.xlsx")