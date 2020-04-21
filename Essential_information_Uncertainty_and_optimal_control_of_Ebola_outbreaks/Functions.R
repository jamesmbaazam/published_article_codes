####This file contains all the functions required to run the 37 models and to simulate the 
####five interventions in the paper titled "Essential information: Uncertainty and optimal 
####control of Ebola outbreaks"

##Gillespie function with Tau-leaping correction 
gillesp <- function(start,ratefun,trans,pars,times=0:50,deltaT) {
  ##Gillespie function with Tau-leaping correction 
  # start = The number of individuals in each compartment at the starting of the epidemic 
  # ratefun = Transition rates between compartments
  # trans = State transition matrix -- gives location and direction of transition
  # pars = List of all the parameters in each model
  # times = The time points at which to store results
  # deltaT = The tau-leap time interval
  
  tfn=function(Y){
    if(tevents>1){
      Z=apply(as.matrix(Y, byrow=T),2,sum)
    }
    if(tevents==1){
      Z=Y
    }
    return(Z)}
  t0 <- times[1]  # set time to starting time              
  ntimes <- length(times)
  X <- start      #	set state to starting state              
  res <- matrix (nrow=length(times),ncol=length(start),
                 dimnames=list(times,names(start)))
  for(ctr in 1:(ntimes-1)) {
    res[ctr,] <- X
    while (t0 < times[ctr+1]) {
      rates <- ratefun(X,pars,t0)
      if(all(rates==0)) break
      totrate <- sum(rates)
      elapsed <- deltaT
      tevents <- rpois(1,totrate*deltaT)
      t0 <- t0 + deltaT # update time
      if(tevents>0){
        which.trans <-
          sample(1:nrow(trans),size=tevents,prob=rates,replace=TRUE) #pick transition     
        while(
          sum((X + tfn(trans[which.trans,]))<0)>0){
          tevents <- max(c(1,rpois(1,totrate*deltaT)))
          which.trans <-
            sample(1:nrow(trans),size=tevents,prob=rates,replace=TRUE)
        }
        X = X + tfn(trans[which.trans,]) # add transition values to current stat
      }
    }
  }
  cbind(times,res)
}

## function to define the transitions between compartments
ratefun <- function(X,pars,time) {
  ## function to define the transitions between compartments
  # pars = List of all the parameters in each model
  # times = The time points at which to store results

  vals <- c(as.list(pars),as.list(X)) # Attach state and parameters as lists
  rates <- with(vals, c(
    #those leave Sc
    ScToEc=(betaIc1*Ic1+betaIc2*Ic2+betaIh1*Ih1+betaIh2*Ih2+betaIccc1*Iccc1+betaIccc2*Iccc2+betaIrh1*Irh1+betaFc*Fc+betaFh*Fh)*Sc/N,
    ScToInccc2=pnccc*(betaIc1*Ic1+betaIc2*Ic2+betaIh1*Ih1+betaIh2*Ih2+betaIccc1*Iccc1+betaIccc2*Iccc2+betaIrh1*Irh1+betaFc*Fc+betaFh*Fh)*Sc/N*alpha*leg*gamma1*ccc1*fcccmax2*Ic1, 
    ScToInh2=pnh*(betaIc1*Ic1+betaIc2*Ic2+betaIh1*Ih1+betaIh2*Ih2+betaIccc1*Iccc1+betaIccc2*Iccc2+betaIrh1*Irh1+betaFc*Fc+betaFh*Fh)*Sc/N*alpha*leg*gamma1*(1-ccc1)*h1*fhmax2*Ic1, 
    ScToSv=rhoV*(Ih1+Ih2+Ihw1)*(Sc/(Sc+Ec+Ic1)), 
    ShwToEhw=(psi*betaIc1*Ic1+psi*betaIc2*Ic2+betaIhwo1*Ihw1+betaIhw1*Ih1+betaIhw2*Ih2+betaIv1*Iv1+betaFhw*Fh+psi*betaFc*Fc)*Shw/N, 
    ShwToSrh=rhoRh*Shw,
    #those leave Srh
    SrhToErh=(betaIc1*Ic1+betaIc2*Ic2+betaIrh1*Irh1+betaFc*Fc)*Srh/N, 
    SrhToShw=rhoH*Srh,
    #those leave Sv
    SvToEv=(betaIhwo1*Ihw1+betaIhw1*Ih1+betaIhw2*Ih2+betaIv1*Iv1+betaFhw*Fh)*Sv/N,   
    SvToSc=rhoRv*Sv,
    #those leave Ec
    EcToIc1=alpha*Ec,    
    EcToEcd=alphaCho*Ec,
    EcToEv=rhoV*(Ih1+Ih2+Ihw1)*(Ec/(Sc+Ec+Ic1)), 
    #those leave Ehw
    EhwToIhw1=alphaHw*Ehw,
    EhwToErh=rhoRh*Ehw,   
    #those leave Erh
    ErhToIrh1=alphaRh*Erh,
    ErhToEhw=rhoH*Erh,   
    #those leave Ev
    EvToIv1=alphaV*Ev,
    EvToEc=rhoRv*Ev,  
    #those leave Ecd    
    EcdToIc1=gammaE1*(1-he1)*Ecd, 
    EcdToIh1=gammaE1*he1*Ecd,
    #those leave Ic1
    Ic1ToIccc1=gamma1*ccc1*fcccmax1*Ic1, 
    Ic1ToIh1=agu*els*gamma1*(1-ccc1)*h1*fhmax1*Ic1,   
    Ic1ToIc2=els*gamma1*(1-ccc1)*(1-h1)*delta1*Ic1, 
    Ic1ToR=gamma1*(1-ccc1)*(1-h1)*(1-delta1)*Ic1,
    Ic1ToIh2=leg*gamma1*(1-ccc1)*h1*fhmax2*Ic1,   
    Ic1ToFc=agu*leg*gamma1*(1-ccc1)*(1-h1)*delta1*Ic1, 
    Ic1ToIccc2=leg*gamma1*ccc1*fcccmax2*Ic1, 
    Ic1ToIv1=rhoV*(Ih1+Ih2+Ihw1)*(Ic1/(Sc+Ec+Ic1)), 
    #those leave Ih1
    Ih1ToIc1=els*gammaH1*phiH1*Ih1, 
    Ih1ToIh2=agu*els*gammaH1*(1-phiH1)*deltaH1*Ih1,
    Ih1ToRh=agu*els*gammaH1*(1-phiH1)*(1-deltaH1)*Ih1, 
    #those leave Ihw1
    Ihw1ToIh2=gammaHw1*hw1*Ihw1,
    Ihw1ToFh=gammaHw1*(1-hw1)*deltaHw1*Ihw1,
    Ihw1ToRh=gammaHw1*(1-hw1)*(1-deltaHw1)*Ihw1, 
    #those leave Irh1
    Irh1ToIhw1=omegaRh*Irh1,
    #those leave Iv1
    Iv1ToIc1=rhoRv*Iv1,
    #those leave Iccc1
    Iccc1ToIc1=gammaCcc1*phiCcc1*Iccc1,
    Iccc1ToFh=mel*gammaCcc1*(1-phiCcc1)*deltaCcc1*Iccc1,
    Iccc1ToIccc2=gammaCcc1*(1-phiCcc1)*deltaCcc1*Iccc1,
    Iccc1ToRccc=gammaCcc1*(1-phiCcc1)*(1-deltaCcc1)*Iccc1,
    #those leave Inccc2
    Inccc2ToEc=gammaNccc*eppsilonNccc*Inccc2,
    Inccc2ToSc=gammaNccc*(1-eppsilonNccc)*Inccc2,
    #those leave Inh2
    Inh2ToEc=gammaNh*eppsilonNh*Inh2,
    Inh2ToSc=gammaNh*(1-eppsilonNh)*Inh2,
    #those leave Ic2
    Ic2ToIccc2=els*gamma2*ccc2*fcccmax2*Ic2, 
    Ic2ToIh2=els*gamma2*(1-ccc2)*h2*fhmax2*Ic2,
    Ic2ToIFc=els*gamma2*(1-ccc2)*(1-h2)*delta2*Ic2,
    Ic2ToR=els*gamma2*(1-ccc2)*(1-h2)*(1-delta2)*Ic2,
    #those leave Ih2
    Ih2ToIc2=els*gammaH2*phiH2*Ih2,
    Ih2ToFh=gammaH2*(1-phiH2)*deltaH2*deltaB2*Ih2,
    Ih2ToR=gammaH2*(1-phiH2)*deltaH2*(1-deltaB2)*Ih2,
    Ih2ToRh=gammaH2*(1-phiH2)*(1-deltaH2)*Ih2,
    #those leave Iccc2
    Iccc2ToIc2=gammaCcc2*phiCcc2*Iccc2,
    Iccc2ToFh=gammaCcc2*(1-phiCcc2)*deltaCcc2*Iccc2,
    Iccc2ToRccc=gammaCcc2*(1-phiCcc2)*(1-deltaCcc2)*Iccc2,
    #those leave Fc;		
    FcToR=gammaFc*Fc,
    #those leave Fh	
    FhToR=gammaFh*Fh,
    #those leave Rh		
    RhToR=gammaRh*Rh,
    #those leave Rccc		
    RcccToR=gammaRccc*Rccc)
  )};


##function to run the 37 models without interventions
runModel <- function(pars.model,modelName,nrun){
# Generic wrapper to run models without interventions
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs

  output <- array(NA, dim=c(length(times),length(save.var),nrun));
  dur=rep(0, nrun); #dur is the duration of epidemic
  pars.t <- pars.model;
  for(r in 1:nrun){
    Gmodel =gillesp(start=c(
      Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
      times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
    output[,,r]=Gmodel[,save.var] 
    dur[r]=sort(output[output[0:(length(times)-1),length(save.var),r]==0,1,r])[1] 
  }  
  save(output, dur, file=paste(modelName,"control.rda",sep='')); 
}


##functions to run the 37 models under the intervention of reducing funeral transmission (abbreviated as "funtra" in the code)
#SEIHFR models
runModel.funtra.SEIHFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHFR models with reduced funeral transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity)) #dur is the duration of epidemic
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaFc"] <- pars.t["betaFc"] * (1- intensity [i]);
    pars.t["betaFh"] <- pars.t["betaFh"] * (1- intensity [i]);
    pars.t["betaFhw"] <- pars.t["betaFhw"] * (1-intensity [i]);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1] 
    }
  }  
  save(output, dur, file=paste(modelName,"funtra.rda",sep=''));
}

#SEIHR models
runModel.funtra.SEIHR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHR models with reduced funeral transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIc1"] <-(1-pCF_F)*pars.t["betaIc1"]+(1- intensity[i])*pCF_F*pars.t["betaIc1"];
    pars.t["betaIc2"] <- (1-pCF_F)*pars.t["betaIc2"]+(1- intensity[i])*pCF_F*pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"funtra.rda",sep='')); 
}

#SEIFR models
runModel.funtra.SEIFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIFR models with reduced funeral transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaFc"] <- pars.t["betaFc"] * (1- intensity[i]);
    pars.t["betaFh"] <- pars.t["betaFh"] * (1- intensity[i]);
    pars.t["betaFhw"] <- pars.t["betaFhw"] * (1-intensity[i]);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"funtra.rda",sep='')); 
}

#SEIR models
runModel.funtra.SEIR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIR models with reduced funeral transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIc1"] <-(1-pF)*pars.t["betaIc1"]+(1- intensity[i])*pF*pars.t["betaIc1"];
    pars.t["betaIc2"] <- (1-pF)*pars.t["betaIc2"]+(1- intensity[i])*pF*pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"funtra.rda",sep=''));
}

##functions to run the 37 models under the intervention of reducing community transmission (abbreviated as "comtra" in the code)
#SEIHFR models
runModel.comtra.SEIHFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHFR models with reduced community transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIc1"] <- pars.t["betaIc1"] * (1- intensity [i]);
    pars.t["betaIc2"] <- pars.t["betaIc2"] * (1- intensity [i]);
    pars.t["betaIrh1"] <- pars.t["betaIrh1"] * (1- intensity [i]);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"comtra.rda",sep='')); 
}

#SEIHR models
runModel.comtra.SEIHR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHR models with reduced community transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIc1"] <-(1-pCF_C)*pars.t["betaIc1"]+(1- intensity [i])*pCF_C*pars.t["betaIc1"];
    pars.t["betaIc2"] <- (1-pCF_C)*pars.t["betaIc2"]+(1- intensity [i])*pCF_C*pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"comtra.rda",sep=''));
}

#SEIFR models
runModel.comtra.SEIFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIFR models with reduced community transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIc1"] <-(1-pCH_C)*pars.t["betaIc1"]+(1- intensity [i])*pCH_C*pars.t["betaIc1"];
    pars.t["betaIc2"] <- (1-pCH_C)*pars.t["betaIc2"]+(1- intensity [i])*pCH_C*pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"comtra.rda",sep=''));
}

#SEIR models
runModel.comtra.SEIR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIR models with reduced community transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIc1"] <-(1-pC)*pars.t["betaIc1"]+(1- intensity [i])*pC*pars.t["betaIc1"];
    pars.t["betaIc2"] <- (1-pC)*pars.t["betaIc2"]+(1- intensity [i])*pC*pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"comtra.rda",sep=''));
}


##functions to run the 37 models under the intervention of reducing mortality ratio (abbreviated as "mor" in the code)
#SEIHFR models
runModel.mor.SEIHFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHFR models with reduced mortality ratio
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["delta1"] <- pars.t["delta1"] * (1-intensity[i]);
    pars.t["delta2"] <- pars.t["delta2"] * (1-intensity[i]);
    pars.t["deltaH1"] <- pars.t["deltaH1"] * (1-intensity[i]);
    pars.t["deltaH2"] <- pars.t["deltaH2"] * (1-intensity[i]);
    pars.t["deltaHw1"] <- pars.t["deltaHw1"] * (1-intensity[i]);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1] 
    }
  }  
  save(output, dur, file=paste(modelName,"mor.rda",sep=''));
}

#SEIHR models
runModel.mor.SEIHR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHR models with reduced mortality ratio
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    wCF_C.new=1/(1+(1-theta.m)*deltaC.m*(1-intensity[i])+theta.m*deltaH.m*(1-intensity[i]))
    wCF_F.new=((1-theta.m)*deltaC.m*(1-intensity[i])+theta.m*deltaH.m*(1-intensity[i]))/(1+(1-theta.m)*deltaC.m*(1-intensity[i])+theta.m*deltaH.m*(1-intensity[i]))
    pars.t["betaIc1"]<-(wCF_C.new/wCF_C)*pCF_C* pars.t["betaIc1"]+(wCF_F.new/wCF_F)*pCF_F* pars.t["betaIc1"];
    pars.t["betaIc2"]<-(wCF_C.new/wCF_C)*pCF_C* pars.t["betaIc2"]+(wCF_F.new/wCF_F)*pCF_F* pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var]
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"mor.rda",sep='')); 
}

#SEIFR models model
runModel.mor.SEIFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIFR models with reduced mortality ratio
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["delta1"] <- pars.t["delta1"] * (1-intensity[i]);
    pars.t["delta2"] <- pars.t["delta2"] * (1-intensity[i]);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"mor.rda",sep='')); 
}

#SEIR model
runModel.mor.SEIR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIR models with reduced mortality ratio
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    wC.new=1/(1+theta.m+(1-theta.m)*deltaC.m*(1-intensity[i])+theta.m*deltaH.m*(1-intensity[i]));
    wH.new=theta.m/(1+theta.m+(1-theta.m)*deltaC.m*(1-intensity[i])+theta.m*deltaH.m*(1-intensity[i]));
    wF.new=((1-theta.m)*deltaC.m*(1-intensity[i])+theta.m*deltaH.m*(1-intensity[i]))/(1+theta.m+(1-theta.m)*deltaC.m*(1-intensity[i])+theta.m*deltaH.m*(1-intensity[i]));
    pars.t["betaIc1"]<-(wC.new/wC)*pC* pars.t["betaIc1"]+(wH.new/wH)*pH* pars.t["betaIc1"]+(wF.new/wF)*pF* pars.t["betaIc1"]
    pars.t["betaIc2"]<-(wC.new/wC)*pC* pars.t["betaIc2"]+(wH.new/wH)*pH* pars.t["betaIc2"]+(wF.new/wF)*pF* pars.t["betaIc2"]
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"mor.rda",sep='')); 
}

##functions to run the 37 models under the intervention of reducing hospital transmission (abbreviated as "hostra" in the code)
#SEIHFR models
runModel.hostra.SEIHFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHFR models with reduced hospital transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIh1"] <- pars.t["betaIh1"] * (1- intensity[i]);
    pars.t["betaIh2"] <- pars.t["betaIh2"] * (1- intensity[i]);
    pars.t["betaIhwo1"] <- pars.t["betaIhwo1"] * (1- intensity[i]);
    pars.t["betaIhw1"] <- pars.t["betaIhw1"] * (1- intensity[i]);
    pars.t["betaIhw2"] <- pars.t["betaIhw2"] * (1- intensity[i]);
    pars.t["betaIv1"] <- pars.t["betaIv1"] * (1- intensity[i]);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1] 
    }
  }  
  save(output, dur, file=paste(modelName,"hostra.rda",sep=''));
}

#SEIHR models
runModel.hostra.SEIHR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHR models with reduced hospital transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model; 
    pars.t["betaIh1"] <- pars.t["betaIh1"] * (1- intensity[i]);
    pars.t["betaIh2"] <- pars.t["betaIh2"] * (1- intensity[i]);
    pars.t["betaIhwo1"] <- pars.t["betaIhwo1"] * (1- intensity[i]);
    pars.t["betaIhw1"] <- pars.t["betaIhw1"] * (1- intensity[i]);
    pars.t["betaIhw2"] <- pars.t["betaIhw2"] * (1- intensity[i]);
    pars.t["betaIv1"] <- pars.t["betaIv1"] * (1- intensity[i]);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"hostra.rda",sep=''));
}

#SEIFR models
runModel.hostra.SEIFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIFR models with reduced hospital transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIc1"] <-(1-pCH_H)*pars.t["betaIc1"]+(1- intensity[i])*pCH_H*pars.t["betaIc1"];
    pars.t["betaIc2"] <- (1-pCH_H)*pars.t["betaIc2"]+(1- intensity[i])*pCH_H*pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"hostra.rda",sep='')); 
}

#SEIR model
runModel.hostra.SEIR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIR models with reduced hospital transmission
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["betaIc1"] <-(1-pH)*pars.t["betaIc1"]+(1- intensity[i])*pH*pars.t["betaIc1"];
    pars.t["betaIc2"] <- (1-pH)*pars.t["betaIc2"]+(1- intensity[i])*pH*pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"hostra.rda",sep=''));
}

  
##functions to run the 37 models under the intervention of increasing hospitalization
#SEIHFR models
runModel.hospitalization.SEIHFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHFR models with reduced increased hospitalization
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["he1"]=ifelse(pars.t["he1"]*(1+intensity[i])<1,pars.t["he1"]*(1+intensity[i]), 1);
    pars.t["h1"]=ifelse(pars.t["h1"]*(1+intensity[i])<1,pars.t["h1"]*(1+intensity[i]), 1); 
    pars.t["h2"]=ifelse(pars.t["h2"]*(1+intensity[i])<1,pars.t["h2"]*(1+intensity[i]), 1); 
    pars.t["hw1"]=ifelse(pars.t["hw1"]*(1+intensity[i])<1,pars.t["hw1"]*(1+intensity[i]), 1);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1] 
    }
  }  
  save(output, dur, file=paste(modelName,"hospitalization.rda",sep='')); 
}

#SEIHR models
runModel.hospitalization.SEIHR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIHR models with reduced increased hospitalization
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    pars.t["he1"]=ifelse(pars.t["he1"]*(1+intensity[i])<1,pars.t["he1"]*(1+intensity[i]), 1); 
    pars.t["h1"]=ifelse(pars.t["h1"]*(1+intensity[i])<1,pars.t["h1"]*(1+intensity[i]), 1); 
    pars.t["h2"]=ifelse(pars.t["h2"]*(1+intensity[i])<1,pars.t["h2"]*(1+intensity[i]), 1); 
    pars.t["hw1"]=ifelse(pars.t["hw1"]*(1+intensity[i])<1,pars.t["hw1"]*(1+intensity[i]), 1);
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"hospitalization.rda",sep='')); 
}

#SEIFR models
runModel.hospitalization.SEIFR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIFR models with reduced increased hospitalization
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    theta.new=ifelse(theta.m*(1+intensity[i])<1,theta.m*(1+intensity[i]),1);
    wCH_C.new=1/(1+theta.new);
    wCH_H.new=theta.new/(1+theta.new);
    pars.t["betaIc1"]<-(wCH_C.new/wCH_C)*pCH_C* pars.t["betaIc1"]+(wCH_H.new/wCH_H)*pCH_H* pars.t["betaIc1"];
    pars.t["betaIc2"]<-(wCH_C.new/wCH_C)*pCH_C* pars.t["betaIc2"]+(wCH_H.new/wCH_H)*pCH_H* pars.t["betaIc2"];
    pars.t["delta1"]<- (1-theta.new)/(1-theta.m)*dC* pars.t["delta1"]+theta.new/theta.m*dH* pars.t["delta1"];
    pars.t["delta2"]<- (1-theta.new)/(1-theta.m)*dC* pars.t["delta2"]+theta.new/theta.m*dH* pars.t["delta2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var]
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"hospitalization.rda",sep=''));
}

#SEIR model
runModel.hospitalization.SEIR <- function(pars.model,modelName,nrun,intensity){
# generic wrapper to run SEIR models with reduced increased hospitalization
# pars.model = a named vector containing model parameters; must match the parameter names required in "modelName"
# modelName = the name of the model, corresponding to the parameter vector in "pars.model"
# nrun = number of stochastic model runs
# intensity = a vector of intensity of intervention

  output <- array(NA, dim=c(length(times),length(save.var),nrun,length(intensity)));
  dur=matrix(NA, nrow=nrun, ncol=length(intensity))
  for( i in 1: length(intensity)){
    pars.t <- pars.model;
    theta.new=ifelse(theta.m*(1+intensity[i])<1,theta.m*(1+intensity[i]),1);
    wC.new=1/(1+theta.new+(1-theta.new)*deltaC.m+theta.new*deltaH.m);
    wH.new=theta.new/(1+theta.new+(1-theta.new)*deltaC.m+theta.new*deltaH.m);
    wF.new=((1-theta.new)*deltaC.m+theta.new*deltaH.m)/(1+theta.new+(1-theta.new)*deltaC.m+theta.new*deltaH.m);
    pars.t["betaIc1"]<-(wC.new/wC)*pC* pars.t["betaIc1"]+(wH.new/wH)*pH* pars.t["betaIc1"]+(wF.new/wF)*pF* pars.t["betaIc1"];
    pars.t["betaIc2"]<-(wC.new/wC)*pC* pars.t["betaIc2"]+(wH.new/wH)*pH* pars.t["betaIc2"]+(wF.new/wF)*pF* pars.t["betaIc2"];
    for(r in 1:nrun){
      Gmodel =gillesp(start=c(
        Sc=9980,Shw=10,Srh=0,Sv=0,Ec=0,Ehw=0,Erh=0,Ev=0,Ecd=0,Ic1=10,Ih1=0,Ihw1=0,Irh1=0,Iv1=0,Iccc1=0,Inccc2=0,Inh2=0,Ic2=0,Ih2=0,Iccc2=0,Fc=0,Fh=0,Rh=0,Rccc=0,R=0,CC=10,CD=0, CP=10),
        times=times,ratefun=ratefun,trans=trans,pars=pars.t,deltaT=0.1);
      cat(date(), " - ", i, "\n")
      output[,,r,i]=Gmodel[,save.var] 
      dur[r,i]=sort(output[output[0:(length(times)-1),length(save.var),r,i]==0,1,r,i])[1]
    }
  }  
  save(output, dur, file=paste(modelName,"hospitalization.rda",sep=''));
}







