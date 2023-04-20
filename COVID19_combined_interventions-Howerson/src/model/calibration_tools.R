## tools to calibrate transmission rate using NGM matix method


# R0 calculations using next-generation matrix method
# function modified from Ottar's book, p. 51
# Istates: vector (character), names of infected classes
# Flist: list (equation), transition rates of new infections
# Vlist: vector (equation), transition rates out of - transition rates into infected classes
# params: vector (double), model parameters
# dfe: vector (double), disease free equilibrium
nextgenR0=function(Istates, Flist, Vlist, params, dfe){
  paras = as.list(c(dfe, params)) 
  k=0
  vl=fl=list(NULL)
  for(i in 1:length(Istates)){
    assign(paste("f", i, sep = "."), lapply(lapply(Flist,deriv, Istates[i]), eval, paras))
    assign(paste("v", i, sep = "."), lapply(lapply(Vlist,deriv, Istates[i]), eval, paras))
    for(j in 1:length(Istates)){
      k=k+1
      fl[[k]]=attr(eval(as.name(paste("f", i, sep=".")))[[j]], "gradient")[1,]
      vl[[k]]=attr(eval(as.name(paste("v", i, sep=".")))[[j]], "gradient")[1,]
    }
  }
  f=matrix(as.numeric(as.matrix(fl)[,1]), ncol=length(Istates))
  v=matrix(as.numeric(as.matrix(vl)[,1]), ncol=length(Istates))
  R0=max(eigen(f%*%solve(v))$values)
  return(R0)
}

### wrapper for nextgenR0
### pass into apply to calc R0 for multiple parameter sets
# in.change.params.names: vector (character), names of change parameters
# in.change.params: vector (double), parameter values to test
# all others: see nextgenR0 description
apply.nextGen = function(in.change.params.names, in.change.params.vals,in.istates, in.flist, in.vlist, in.params,in.dfe, in.sim){
  new.params = in.params
  for(i in 1:length(in.change.params.vals)){
    name = as.character(in.change.params.names[i])
    val = in.change.params.vals[i]
    new.params[name] = val
  }
  tmp = nextgenR0(Istates=in.istates, Flist=in.flist, Vlist=in.vlist, params=new.params, dfe=in.dfe)
  if(in.sim %% 10000 == 0) {print(in.sim)}
  ret = in.change.params.vals
  names(ret) = in.change.params.names
  ret = as.data.frame(t(ret))
  return(data.frame(ret,R0 = tmp))
}

# implement covid model (with testing) setup
# for NGM calculations
mod_config = function(){
  istates = c("E","I1","I2","I3","I4","A","TI1","TI2","TI3","TI4","TA","H") #,"TI2"
  ## new infections
  flist = c(dEdt = quote((1-d)*beta*(I1+I2+I3+I4+rho*A)*S/(S+E+I1+I2+I3+I4+A+R+TR)), # tranmission parameters
            dI1dt = quote(0), dI2dt = quote(0),
            dI3dt = quote(0), dI4dt = quote(0), dAdt = quote(0),
            dTI1dt = quote(0), dTI2dt = quote(0),
            dTI3dt = quote(0), dTI4dt = quote(0), dTAdt = quote(0), dHdt = quote(0))  #,dTI2dt = quote(0)
  ## all lost transfers (out of infected compartments)
  Vm1 = quote(theta*E)                               # dE/dt = to presymptomatic and asymptomatic
  Vm2 = quote((gamma_I1I2+Ta)*I1)                    # dI1/dt = testing I + reported symptoms
  Vm3 = quote((gamma_I2R +Ts)*I2)                    # dI2/dt = testing A + recovery 
  Vm4 = quote((gamma_I1I2+Ta)*I3)                    # dI3/dt = testing I3 + to I4
  Vm5 = quote((gamma_I4H+Ts)*I4)                     # dI4/dt = to hospital
  Vm6 = quote((gamma_AR+Ta)*A)                       # dA/dt = testing P + to symptomatic
  Vm7 = quote(gamma_I1I2*TI1)                        # dTI1/dt = to TI2
  Vm8 = quote(gamma_I2R*TI2)                         # dTI2/dt = to TR
  Vm9 = quote(gamma_I1I2*TI3)                        # dTI3/dt = to TI4
  Vm10 = quote(gamma_I4H*TI4)                        # dTI4/dt = to H
  Vm11 = quote(gamma_AR*TA)                          # dTA/dt = to TR
  Vm12 = quote((gamma_HTR+alpha)*H)                  # dH/dt = to TR + to D
  ## all gains into each infected compartment (that are not new infections -- transfers among infected classes)
  Vp1 = quote(0)                                     # dE/dt = none - all new infections
  Vp2 = quote(((1-p)*(1-q)*theta*E))                 # dI1/dt = from E
  Vp3 = quote(gamma_I1I2*I1)                         # dI2/dt = from I1
  Vp4 = quote((1-p)*q*theta*E)                       # dI3/dt = from E
  Vp5 = quote(gamma_I1I2*I3)                         # dI4/dt = fromm I3
  Vp6 = quote(p*theta*E)                             # dA/dt = from E
  Vp7 = quote(Ta*I1)                                 # dTI1/dt = from testing
  Vp8 = quote((Ts*I2)+gamma_I1I2*TI1)                # dTI2/dt = from testing + from TI1
  Vp9 = quote(Ta*I3)                                 # dTI3/dt = from testing
  Vp10 = quote(Ts*I4 + gamma_I1I2*TI3)               # dTI4/dt = from testing + from TI3
  Vp11 = quote(Ta*A)                                 # dTA/dt = from testing
  Vp12 = quote(gamma_I4H*(TI4+I4))                   # dH/dt = from I4 + from TI4
  ## Vm - Vp
  V1 = substitute(a - b, list(a = Vm1, b = Vp1))
  V2 = substitute(a - b, list(a = Vm2, b = Vp2))
  V3 = substitute(a - b, list(a = Vm3, b = Vp3))
  V4 = substitute(a - b, list(a = Vm4, b = Vp4))
  V5 = substitute(a - b, list(a = Vm5, b = Vp5))
  V6 = substitute(a - b, list(a = Vm6, b = Vp6))
  V7 = substitute(a - b, list(a = Vm7, b = Vp7))
  V8 = substitute(a - b, list(a = Vm8, b = Vp8))
  V9 = substitute(a - b, list(a = Vm9, b = Vp9))
  V10 = substitute(a - b, list(a = Vm10, b = Vp10))
  V11 = substitute(a - b, list(a = Vm11, b = Vp11))
  V12 = substitute(a - b, list(a = Vm12, b = Vp12))
  ## make vlist
  vlist = c(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12)
  ## disease free equilibrium
  df = list(S=1, E=0, I1=0, I2=0,I3=0,I4=0,A=0,TI1=0,TI2 = 0,TI3=0,TI4=0,TA=0,H=0,R=0,TR=0)
  return(list(istates = istates, flist = flist, vlist = vlist, dfe = df))
}