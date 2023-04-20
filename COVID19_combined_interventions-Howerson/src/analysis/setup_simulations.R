
###### Read command line args ##########
args <- commandArgs(trailingOnly = TRUE)

input_params <- args[1]
input_calibration_tools <- args[2]
output_file <- args[3]


# input_params <- "output/sims/parameters.rda" 
# input_calibration_tools <- "src/model/calibration_tools.R"
# output_file <- "output/sims/sims_setup.rda"


source(input_calibration_tools)
params = readRDS(input_params)
list2env(params, globalenv())

output = list()

# setting up for following sim types:
#    main = base intervention simulations with best estimate parameters
#    sensitivity_IC = simulations with varying initial conditions
#    sensitivity_asymp = simulations with varying assumptions about % asymptomatic infections
#    sensitivity_asymp_rho = simulations with varying assumptions about rel infectiousness of asymptomatic infections


## main
tst.intervention = expand.grid(d = test_d,
                               tot.tests = test_tottests,
                               sens = test_sens,
                               Ds = test_Ds)
tst.intervention$Da = tst.intervention$Ds
tst.intervention$Ts = 1/tst.intervention$Ds
tst.intervention$Ta = 1/tst.intervention$Da
tst.intervention$sim = 1:nrow(tst.intervention)

tst.asymp = data.frame(p = params["p"])
tst.asymp$current.sim = 1:nrow(tst.asymp)

# save
output[["main"]] = list(int = tst.intervention, asymp = tst.asymp)
# intervention$main = tst.intervention
# asymp$main = tst.asymp


## sensitivity_IC
tst.intervention = expand.grid(d = test_d,
                               tot.tests = test_tottests,
                               sens = 1,
                               Ds = test_Ds)
tst.intervention$Da = tst.intervention$Ds
tst.intervention$Ts = 1/tst.intervention$Ds
tst.intervention$Ta = 1/tst.intervention$Da
tst.intervention$sim = 1:nrow(tst.intervention)

tst.asymp = data.frame(p = params["p"])
tst.asymp$current.sim = 1:nrow(tst.asymp)

# save
output[["sensitivity_IC"]] = list(int = tst.intervention, asymp = tst.asymp)
# intervention$sensitivity_IC = tst.intervention
# asymp$sensitivity_IC = tst.asymp



## sensitivity_asymp
tst.intervention = expand.grid(d = test_d,
                               tot.tests = test_tottests,
                               sens = 1,
                               Ds = test_Ds,
                               p = c(0.2,0.6))
tst.intervention$Da = tst.intervention$Ds
tst.intervention$Ts = 1/tst.intervention$Ds
tst.intervention$Ta = 1/tst.intervention$Da
tst.intervention$sim = 1:nrow(tst.intervention)


asymp = data.frame(p = c(0.2,0.6), beta = NA)

params.vary = params
mod = mod_config()
for(j in 1:nrow(asymp)){
  #browser()
  params.vary["p"] = asymp$p[j]
  rslts = data.frame(beta = seq(0.4,0.6, by = 0.001))
  rslts = apply(rslts, 1, apply.nextGen, in.change.params.names = c("beta"),
                in.istates = mod$istates, in.flist = mod$flist, in.vlist = mod$vlist,
                in.params = params.vary, in.dfe = mod$df, in.sim = 1)
  rslts = do.call(rbind, rslts)
  est.beta = data.frame(R0 = c(2.5), beta = NA)
  for(i in est.beta$R0){
    est.beta$beta[which(est.beta$R0 == i)] = rslts$beta[which(abs(rslts$R0-i) == min(abs(rslts$R0-i)))]
  }
asymp$beta[j] = est.beta$beta[est.beta$R0 == 2.5]
}

tst.asymp = data.frame(p = c(0.2,0.6))
tst.asymp = merge(tst.asymp, asymp)
tst.asymp$current.sim = 1:nrow(tst.asymp)

# save
output[["sensitivity_asymp"]] = list(int = tst.intervention, asymp = tst.asymp)
# intervention$sensitivity_asymp = tst.intervention
# asymp$sensitivity_asymp = tst.asymp


## sensitivity_asymp_rho
tst.intervention = expand.grid(d = test_d,
                               tot.tests = test_tottests,
                               sens = 1,
                               Ds = test_Ds,
                               rho = 0.5)
tst.intervention$Da = tst.intervention$Ds
tst.intervention$Ts = 1/tst.intervention$Ds
tst.intervention$Ta = 1/tst.intervention$Da
tst.intervention$sim = 1:nrow(tst.intervention)

params.vary = params
params.vary["rho"] = 0.5
mod = mod_config()
rslts = data.frame(beta = seq(0.4,0.8, by = 0.001))
rslts = apply(rslts, 1, apply.nextGen, in.change.params.names = c("beta"),
                in.istates = mod$istates, in.flist = mod$flist, in.vlist = mod$vlist,
                in.params = params.vary, in.dfe = mod$df, in.sim = 1)
rslts = do.call(rbind, rslts)

tst.asymp = data.frame(p = params["p"], beta = rslts$beta[which(abs(rslts$R0-2.5) == min(abs(rslts$R0-2.5)))])
tst.asymp$current.sim = 1:nrow(tst.asymp)

# save
output[["sensitivity_asymp_rho"]] = list(int = tst.intervention, asymp = tst.asymp)
# intervention$sensitivity_asymp_rho = tst.intervention
# asymp$sensitivity_asymp_rho = tst.asymp


saveRDS(output, output_file)