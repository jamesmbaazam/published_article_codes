### step fn to implement SEIR model
### for implemeting limited numbers of tests available per day
# sims: integer, number of stochastic replicates
# delta.t: double, length of time step
# S,E,I1,...N,D: vector, number of inidivudauls in class i at time t
# all else: double, parameter values
fixedtest_step <- function (sims,delta.t,S,E,I1,I2,I3,I4,A,R,H,
                            WI1, WI2, WI3, WI4, WA,
                            TI1,TI2,TI3,TI4,TA,N,D, rho,
                            d,beta,theta,p,q,gamma_I1I2,gamma_I2R, gamma_HD,
                            gamma_AR,gamma_I4H, gamma_HTR, 
                            Ts, Ta, tot.test, sens,spec,c, ...) {
  # adapted from Aaron King's code
  ### transitions 
  # unisolated classes
  dN_SE <- rbinom(n=sims,size=S,prob=1-exp(-1*(1-d)*beta*(I2+I4+WI2+WI4+rho*(I1+I3+A+WI1+WI3+WA))/(N-TI1-TI2-TI3-TI4-TA-H-D)*delta.t))
  dN_EI <- rbinom(n=sims,size=E,prob=1-exp(-theta*delta.t))
  dN_EI1 <- rbinom(n=sims,size=dN_EI,prob=(1-p)) # number of symptomatics
  dN_EA <- dN_EI - dN_EI1
  dN_EI3 <- rbinom(n = sims, size = dN_EI1, prob = q) # number of hosp. bound cases
  dN_EI1 <- dN_EI1 - dN_EI3
  dN_I1I2 <- rbinom(n=sims,size=I1,prob=1-exp(-gamma_I1I2*delta.t))
  dN_I3I4 <- rbinom(n=sims,size=I3,prob=1-exp(-gamma_I1I2*delta.t)) #assuming I1->I2 = I3->I4
  dN_I2R <- rbinom(n=sims,size=I2,prob=1-exp(-gamma_I2R*delta.t))
  dN_AR <- rbinom(n = sims,size =A,prob = 1-exp(-gamma_AR*delta.t))
  dN_I4H <- rbinom(n = sims, size = I4, prob = 1-exp(-gamma_I4H*delta.t))
  dN_HR <- rbinom(n = sims, size = H, prob = 1-exp(-gamma_HTR*delta.t))
  # waiting for test results 
  dN_WI1WI2 <- rbinom(n=sims,size=WI1,prob=1-exp(-gamma_I1I2*delta.t))
  dN_WI3WI4 <- rbinom(n=sims,size=WI3,prob=1-exp(-gamma_I1I2*delta.t))
  dN_WI4H <- rbinom(n=sims,size=WI4,prob=1-exp(-gamma_I4H*delta.t))
  dN_WI2R <- rbinom(n=sims,size=WI2,prob=1-exp(-gamma_I2R*delta.t))
  dN_WAR <- rbinom(n = sims,size = WA,prob = 1-exp(-gamma_AR*delta.t))
  # tested classes
  dN_TI1TI2 <- rbinom(n=sims,size=TI1,prob=1-exp(-gamma_I1I2*delta.t))
  dN_TI3TI4 <- rbinom(n=sims,size=TI3,prob=1-exp(-gamma_I1I2*delta.t))
  dN_TI4H <- rbinom(n=sims,size=TI4,prob=1-exp(-gamma_I4H*delta.t))
  dN_TI2R <- rbinom(n=sims,size=TI2,prob=1-exp(-gamma_I2R*delta.t))
  dN_TAR <- rbinom(n = sims,size = TA,prob = 1-exp(-gamma_AR*delta.t))
  ### update states
  # unisolated classes
  S <- S - dN_SE
  E <- E + dN_SE - dN_EI1 - dN_EA - dN_EI3
  I1 <- I1 + dN_EI1 - dN_I1I2
  I2 <- I2 + dN_I1I2 - dN_I2R
  I3 <- I3 + dN_EI3 - dN_I3I4
  I4 <- I4 + dN_I3I4 - dN_I4H
  A <- A + dN_EA - dN_AR
  R <- R + dN_I2R + dN_AR + dN_TI2R + dN_TAR + dN_HR + dN_WAR + dN_WI2R
  H <- H + dN_I4H + dN_TI4H - dN_HR + dN_WI4H
  #D <- D + dN_HD
  # waiting for test results 
  WI1 <- WI1 - dN_WI1WI2
  WI2 <- WI2 + dN_WI1WI2 - dN_WI2R
  WI3 <- WI3 - dN_WI3WI4
  WI4 <- WI4 + dN_WI3WI4 - dN_WI4H
  WA <-  WA  - dN_WAR
  # isolated classes
  TI1 <- TI1 - dN_TI1TI2
  TI2 <- TI2 + dN_TI1TI2 - dN_TI2R
  TI3 <- TI3 - dN_TI3TI4
  TI4 <- TI4 + dN_TI3TI4 - dN_TI4H
  TA <- TA - dN_TAR
  # counters
  new_cases = dN_SE
  new_hosp = dN_I4H + dN_TI4H + dN_WI4H
  new_symp = dN_I1I2 + dN_I3I4
  ### TEST RESULTS
  # deaths 
  dN_HD <-  rbinom(n = sims, size = H, prob = 1-exp(-gamma_HD*delta.t))
  #browser()
  H <- H - dN_HD
  D <- D + dN_HD
  ### TESTING
  test.symp = matrix(0, sims, 2)
  test.asymp = matrix(0,sims, 3)
  test.symp.numpos = rep(0,sims)
  extras = c*(1+p)*N # assume some % of the population will present every day with symptoms and not be infected with COVID, this will increase with increased reporting
  num.symp.test = pmax_fast(0,pmin_fast(round(I2+I4+extras),tot.test))
  num.asymp.test = tot.test - num.symp.test
  ### symptomatic tests
  if(any(num.symp.test>0)){
    symp.to.test = which((I2+I4)>0)
    test.symp.numpos[symp.to.test] = rbinom(n = length(symp.to.test), size = num.symp.test, prob = (sens*(I2+I4)/(I2+I4+extras))[symp.to.test])
    test.symp = mapply(symp_sample,I2,I4,extras, test.symp.numpos, USE.NAMES = TRUE)
    test.symp = t(test.symp)
  }
  ### asymptoamtic tests
  if(any(num.asymp.test>0)){
    test.asymp.numpos = rbinom(n = sims, size = num.asymp.test, prob = sens*(I1+I3+A)/(S+E+I1+A+R+I3))
    test.asymp = mapply(asymp_sample, S = S, E = E, I1 = I1,I3 = I3,A = A,R = R, test.asymp.numpos, USE.NAMES = TRUE)
    test.asymp = t(test.asymp)
  }
  I1 <- I1 - test.asymp[,1]
  I2 <- I2 - test.symp[,1]
  I3 <- I3 - test.asymp[,2]
  I4 <- I4 - test.symp[,2]
  A <- A - test.asymp[,3]  
  # isolated, testing classes
  WI1 <- WI1 + test.asymp[,1]
  WI2 <- WI2 + test.symp[,1]
  WI3 <- WI3 + test.asymp[,2]
  WI4 <- WI4 + test.symp[,2]
  WA <- WA + test.asymp[,3]
  # transition to isolation
  dN_WI1TI1 <- rbinom(n=sims,size=WI1,prob=1-exp(-Ta*delta.t))
  dN_WI2TI2 <- rbinom(n=sims,size=WI2,prob=1-exp(-Ts*delta.t))
  dN_WI3TI3 <- rbinom(n=sims,size=WI3,prob=1-exp(-Ta*delta.t))
  dN_WI4TI4 <- rbinom(n=sims,size=WI4,prob=1-exp(-Ts*delta.t))
  dN_WATA <- rbinom(n=sims,size=WA,prob=1-exp(-Ta*delta.t))
  # waiting for test results 
  WI1 <- WI1 -  dN_WI1TI1
  WI2 <- WI2 - dN_WI2TI2
  WI3 <- WI3 - dN_WI3TI3
  WI4 <- WI4 - dN_WI4TI4
  WA <- WA - dN_WATA
  # isolated classes
  TI1 <- TI1 + dN_WI1TI1
  TI2 <- TI2 + dN_WI2TI2
  TI3 <- TI3 + dN_WI3TI3
  TI4 <- TI4 + dN_WI4TI4
  TA <- TA + dN_WATA
  # counters
  N = S+E+I1+I2+I3+I4+A+R+H+TI1+TI2+TI3+TI4+TA+WI1+WI2+WI3+WI4+WA+D
  ret = cbind(S,E,I1,I2,I3,I4,A,R,H,WI1,WI2,WI3,WI4,WA,
              TI1,TI2,TI3,TI4,TA,D,N,new_cases,new_hosp, num.symp.test, new_sympiso = (dN_WI2TI2+dN_WI4TI4), new_asympiso = (dN_WI1TI1+dN_WI3TI3+dN_WATA)) 
  if(any(ret<0)){browser();print(ts);ret = update.negs(ret)}
  return(ret)
}