# ----- RUN MODELS ----- #

## function for asymptomatic sampling
## given test.asymp.numpos positive asymptoatmic tests, how many of those are from I1,I3, and A classes
# pass in single values from each class (use apply for multiple sims)
asymp_sample = function(S,E,I1,I3,A,R, test.asymp.numpos){
  if(test.asymp.numpos == 0){return(c(0,0,0))}
  if(test.asymp.numpos >= (I1+ I3 + A)){return(c(I1,I3,A))}
  m <- c(I1,I3,A)
  # expand the item counts in to a single vector with i repeated m[i] times
  probs <- c(rep(m[1],m[1]),rep(m[2],m[2]),rep(m[3],m[3]))/(I1+I3+A+S+E+R)
  s = sample_int_ccrank(sum(m), test.asymp.numpos, probs)
  return(c(length(s[s<=m[1]]),length(s[s<=(m[2]+m[1]) & s>m[1]]),length(s[s>(m[2]+m[1])])))
}

## function for symptomatic sampling
## given symp.test available symptomatic tests, how many of those are positive from I2 and I4 classes
# pass in single values from each class (use apply for multiple sims)
# extras represents other individuals that report symptoms, but don't have covid
symp_sample = function(I2,I4,extras, test.symp.numpos){
  if(test.symp.numpos == 0){return(c(0,0))}
  if(test.symp.numpos >= (I2 + I4)){return(c(I2,I4))}
  m <- c(I2,I4)
  # expand the item counts in to a single vector with i repeated m[i] times
  probs <- c(rep(m[1],m[1]),rep(m[2],m[2]))/(I2+I4+extras)
  s = sample_int_ccrank(sum(m), test.symp.numpos, probs)
  return(c(length(s[s<=m[1]]),length(s[s<=(m[2]+m[1]) & s>m[1]])))
}


### create container for simulation output, add initial conditions
# inits: vector (integer), initial values for each class
# sims: number of stochastic replicates
fixedtest_setup_IC = function(inits,sims){
  with(as.list(inits),{
    NC = ifelse("NC_Diff" %in% names(inits), NC_Diff, new_cases)
    IC = cbind(S = rep(S, sims), E = rep(E,sims),
               I1 = rep(I1, sims), I2 = rep(I2,sims),
               I3 = rep(I3,sims), I4 = rep(I4,sims),
               A = rep(A, sims), R = rep(R,sims), 
               H = rep(H, sims),WI1 = rep(0,sims),
               WI2 = rep(0, sims), WI3 = rep(0,sims),
               WI4 = rep(0,sims), WA = rep(0,sims),
               TI1 = rep(TI1,sims),TI2 = rep(TI2, sims), 
               TI3 = rep(TI3,sims),TI4 = rep(TI4,sims), 
               TA = rep(TA,sims), 
               D = rep(0,sims),
               N = rep(S+E+I1+I2+I3+I4+A+R+H+TI1+TI2+TI3+TI4+TA, sims), 
               new_cases = rep(NC, sims), new_hosp = rep(0, sims), symp_tests = rep(0,sims), new_sympiso = rep(0, sims), new_asympiso = rep(0, sims))#,
    # new_symp = rep(0, sims),
    # true_pos = rep(0,sims), true_neg = rep(0,sims),
    # false_pos = rep(0, sims), false_neg = rep(0,sims))
    return(IC)
  })
}

### implement fixedtest_step over time
# IC: vector, initial conditions for each state variable
# params: vector, parameter values
# times: integer, number of time steps to run for
# change.vars: any parameters to update
fixedtest_run_sims = function(IC, params, times, change.vars){
  with(as.list(c(change.vars,params, IC)),{
    class = colnames(IC)
    sims = nrow(IC)
    ret = array(integer(),c(sims,length(class),times),dimnames = list(NULL, class, NULL))
    # ret dimensions - 1: sims, 2: classes, 3: time
    ret[,,1] = IC
    # run first phase of intervention (distanced, no testing)
    for(ts in 2:times){
      ret[,,ts] <- fixedtest_step(sims,delta.t =1, S = ret[,"S",ts-1], E = ret[,"E",ts-1], 
                                  I1 = ret[,"I1",ts-1],I2 = ret[,"I2",ts-1],A = ret[,"A",ts-1],
                                  I3 = ret[,"I3",ts-1],I4 = ret[,"I4",ts-1], R = ret[,"R",ts-1],
                                  H = ret[,"H",ts-1],WI1 = ret[,"WI1",ts-1],WI2 = ret[,"WI2",ts-1], 
                                  WI3 = ret[,"WI3",ts-1],WI4 = ret[,"WI4",ts-1],WA = ret[,"WA",ts-1],
                                  TI1 = ret[,"TI1",ts-1],TI2 = ret[,"TI2",ts-1], #,D = ret[,"D",ts-1]
                                  TI3 = ret[,"TI3",ts-1],TI4 = ret[,"TI4",ts-1],TA = ret[,"TA",ts-1],
                                  N = ret[,"N",ts-1],D = ret[,"D",ts-1],
                                  d= d,beta=beta,theta=theta,p=p,q=q,c=c,rho = rho,
                                  gamma_I1I2=gamma_I1I2,gamma_I2R=gamma_I2R,#alpha = alpha,
                                  gamma_AR=gamma_AR,gamma_I4H=gamma_I2H, gamma_HTR=gamma_HTR, gamma_HD = gamma_HD,
                                  Ts=Ts,Ta = Ta,tot.test= tot.tests,sens=sens,spec=spec)
    }
    return(list(ret))
  })
}

### implement fixedtest_step over time, return objective quantiles only (for speed)
# IC: vector, initial conditions for each state variable
# params: vector, parameter values
# times: integer, number of time steps to run for
# change.vars: any parameters to update
fixedtest_run_sims_return = function(IC, params, times, change.vars){
  with(as.list(c(change.vars,params, IC)),{
    #if(!any(is.na(change.vars))){if(sim %% 10000 == 0){print(paste("sim:",sim,"time:",Sys.time()-start_time))}}
    class = colnames(IC)
    sims = nrow(IC)
    ret = array(integer(),c(sims,length(class),times),dimnames = list(NULL, class, NULL))
    # ret dimensions - 1: sims, 2: classes, 3: time
    ret[,,1] = IC
    # run first phase of intervention (distanced, no testing)
    for(ts in 2:times){
      d.tmp = ifelse(ts == 2, 0.5, d)
      #    d.tmp = d
      ret[,,ts] <- fixedtest_step(sims,delta.t =1, S = ret[,"S",ts-1], E = ret[,"E",ts-1], 
                                  I1 = ret[,"I1",ts-1],I2 = ret[,"I2",ts-1],A = ret[,"A",ts-1],
                                  I3 = ret[,"I3",ts-1],I4 = ret[,"I4",ts-1], R = ret[,"R",ts-1],
                                  H = ret[,"H",ts-1],WI1 = ret[,"WI1",ts-1],WI2 = ret[,"WI2",ts-1], 
                                  WI3 = ret[,"WI3",ts-1],WI4 = ret[,"WI4",ts-1],WA = ret[,"WA",ts-1],
                                  TI1 = ret[,"TI1",ts-1],TI2 = ret[,"TI2",ts-1], #,D = ret[,"D",ts-1]
                                  TI3 = ret[,"TI3",ts-1],TI4 = ret[,"TI4",ts-1],TA = ret[,"TA",ts-1],
                                  N = ret[,"N",ts-1],D = ret[,"D",ts-1],
                                  d= d.tmp,beta=beta,theta=theta,p=p,q=q,c= c,rho =rho,
                                  gamma_I1I2=gamma_I1I2,gamma_I2R=gamma_I2R,#alpha = alpha,
                                  gamma_AR=gamma_AR,gamma_I4H=gamma_I2H, gamma_HTR=gamma_HTR, gamma_HD = gamma_HD,
                                  Ts=Ts,Ta = Ta, tot.test= tot.tests,sens=sens,spec=spec)
    }
    return(list(sim_return_quant(ret)))
  })
}

### shell fn to implement fixedtest_run_sims over multiple parameter values
### pass into apply parameters that need to be same for initial lockdown runs too
### (e.g. p needs to be consistent for lockdown and post-lockdown, d does not)
# change.vars: matrix, parameter combinations to be tested in post-lockdown
# params: vector, overall params (should be parameters for lockdown period)
# times.first: length of lockdown period
# times.second: length of post-lockdown period
# inits: vector, initial conditions for each state variable
# sims: integer, number of stochastic replicates
fixedtest_apply = function(change.vars,params,IC,times.first,times.second,inits,sims){ #return
  with(as.list(c(inits,params)),{
    #if(current.sim%%50==1){print(current.sim)}
    print(paste("p:",p))
    new.params = params
    for(i in names(inits)[names(inits)%in% names(params)]){
      new.params[i] = inits[i] # update parameters/ICs with new values (stored in inits)
      change.vars = change.vars %>% filter(get(i) == new.params[i]) # subset change.vars to only those relevant to current parameter set
    }
    delta.t = 1
    #setup IC
    #IC = fixedtest_setup_IC(inits, sims)
    # run first phase of intervention (distanced, no testing)
    #ret.first = fixedtest_run_sims(IC, new.params, times.first, change.vars = NA)
    # set IC for second phase as end of first phase
    #IC = ret.first[[1]][,,times.first]
    # run second phase over multiple scenarios
    # out = apply(change.vars,1,fixedtest_run_sims, IC = IC,times = times.second, params=new.params) #ret.second =
    # out = do.call(rbind, out)
    # out = lapply(1:length(out), function(x){cbind(change.vars[x,], melt(out[[x]][,"new_cases",]))})
    # out = do.call(rbind, out)
    out = apply(change.vars,1,fixedtest_run_sims_return, IC = IC,times = times.second, params=new.params) #ret.second =
    out = do.call(rbind, out)
    out = do.call(rbind, out)
    out = cbind(change.vars[rep(seq_len(nrow(change.vars)), each = 7), ], as.data.frame(out))
    return(list(out))
  })
}

# ----- SIMULATION UTILITIES ----- #

update.negs = function(ret){
  browser()
  neg = which(ret<0 | is.na(ret)) # identify negative cells
  neg = rbind(neg,ret[neg]) # amount they are negative
  print(neg[2,])
  neg = rbind(neg, ifelse(neg[1,] %% nrow(ret)==0, nrow(ret),neg[1,] %% nrow(ret))) # row of each neg 
  ret[neg[1,]] = ret[neg[1,]] - neg[2,] # remove negative from inf. col
  neg = rbind(neg, ret[neg[3,],"R"]) # number of recovered
  ret[neg[3,],"R"] = ret[neg[3,],"R"] + neg[2,] # take out of recovered class
  return(ret)
}

### function to calculate objectives from time series
### objectives: total cases, max hosp., days abve hosp. thresh, deaths
### pass in single time/vars matrix, then apply this function across all sims in the array (so returning outcomes for all sims)
# dat: matrix (double), simulation output (number of individuals in each class over time)
# H.thresh: integer, hospital threshold
calc.obj = function(dat, H.thresh){
  #browser()
  ret = vector(length = 5, mode = "integer")
  # tot cases
  ret[1] = sum(dat["new_cases",])
  # peak hosp
  ret[2] = max(dat["H",])
  # days above thresh
  ret[3] = length(which(dat["H",]>H.thresh))
  # total hosp
  #browser()
  ret[4] = sum(dat["new_hosp",])
  # deaths
  ret[5] = dat["D",ncol(dat)]
  # total symp isos
  ret[6] = sum(dat["new_sympiso",])
  # total asymp isos
  ret[7] = sum(dat["new_asympiso",])
  names(ret) = c("tot.cases", "peak.hosp", "days.above.thresh", "tot.hosp","death", "tot.sympsio", "tot.asympiso")
  return(ret)
}

### return 5,50,95 quantiles of objective outcomes
# dat: matrix, model output
sim_return_quant = function(dat){
  dat = apply(dat, 1, calc.obj, H.thresh = H.thresh)
  quant = matrix(integer(), 3,5)
  quant.1 = apply(dat, 1, median)
  quant.2 = apply(dat,1,quantile, 0.05)
  quant.3 = apply(dat,1,quantile, 0.95)
  quant = rbind(quant.1,quant.2,quant.3)
  return(t(quant))
}

### find elementwise min/max, updated for speed
# from https://r.789695.n4.nabble.com/pmin-pmax-slower-than-necessary-for-common-cases-td909143.html
# k,x: vectors
pmin_fast = function(k,x) {(x+k - abs(x-k))/2}
pmax_fast = function(k,x) {(x+k + abs(x-k))/2}
