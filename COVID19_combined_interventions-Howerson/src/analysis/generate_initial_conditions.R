## Preamble
suppressWarnings({
  suppressMessages({
    require(knitr)
    require(ggplot2)
    require(tidyr)
    require(gridExtra)
    require(ggplot2)
    require(gridExtra)
    require(RColorBrewer)
    require(plyr)
    require(dplyr)
    require(deSolve)
    require(pmultinom)
    require(abind)
    require(cowplot)
    require(reshape)
    require(reshape2)
    require(metR)
    require(grid)
    require(scales)
    require(viridis)
    require(wrswoR)
  })})



###### Read command line args ##########
args <- commandArgs(trailingOnly = TRUE)

input_model_file <- args[1]
input_model_tools <- args[2]
input_parameters <- args[3]
output_file <- args[4]

# input_model_file <- "src/model/model_immediate_symp_isolation.R"
# input_model_tools <- "src/model/model_tools.R"
# input_parameters <- "output/sims/parameters.rda"
# output_initial_conditions_full <- "output/sims/initial_sims_full.rda"
# output_initial_conditions_low <- "output/sims/initial_conditions_low.csv"
# output_initial_conditions_mid <- "output/sims/initial_conditions_mid.csv"
# output_initial_conditions_high <- "output/sims/initial_conditions_high.csv"



##### generate initial conditions #####

source(input_model_file)
source(input_model_tools)

params = readRDS(input_parameters)
list2env(params, globalenv())


# setup initial conditions
# all stochastic runs
sims = 1000


# initial conditions
start = c(S = 9998, E = 0, I1 = 1, I2 = 0, I3 = 0, I4 = 0,
          A = 1, R = 0, H = 0, D = 0,
          WI1 = 0, WI2 = 0, WI3 = 0, WI4 = 0, WA = 0,
          TI1 = 0, TI2 = 0, TI3 = 0,TI4 = 0, TA = 0, new_cases = 2)

# phase 1: unbounded growth
inits = fixedtest_run_sims(IC = fixedtest_setup_IC(start, sims), params, t_unbounded, change.vars = NA)
inits = inits[[1]]
# phase 2: transmission slowed by NPIs
inits_2 = fixedtest_run_sims(IC = inits[,,t_unbounded], params, t_distanced, change.vars = c( d = d_distanced))
inits_2 = inits_2[[1]]
inits = abind(inits, inits_2)
# remove those that faded out
fadeouts = which(apply(inits[,c("I1", "I2","I3","I4","A"),t_unbounded + t_distanced], 1,sum)==0)
n_fadeouts = length(fadeouts)
if(n_fadeouts>0){inits = inits[-fadeouts, ,]}
# recalculate 
while(n_fadeouts>0){
  #browser()
  inits_tmp = fixedtest_run_sims(IC = fixedtest_setup_IC(start, n_fadeouts), params, t_unbounded, change.vars = NA)
  inits_tmp = inits_tmp[[1]]
  new_IC = if(dim(inits)[1] == sims- 1){new_IC = t(as.matrix(inits_tmp[,,t_unbounded]))}
  else{new_IC = inits_tmp[,,t_unbounded]}
  inits_tmp_2 = fixedtest_run_sims(IC = new_IC, params, t_distanced, change.vars = c(d = d_distanced))
  inits_tmp_2 = inits_tmp_2[[1]]
  inits_tmp = abind(inits_tmp, inits_tmp_2)
  if(dim(inits)[1] == sims - 1){
    if(sum(inits_tmp[,c("I1", "I2","I3","I4","A"),t_unbounded + t_distanced])==0){n_fadeouts = 1}
    else{n_fadeouts = 0}}
  if(dim(inits)[1] != sims - 1){
    fadeouts_tmp = which(apply(inits_tmp[,c("I1", "I2","I3","I4","A"),t_unbounded + t_distanced], 1,sum)==0);
    n_fadeouts = length(fadeouts_tmp)}
  if(n_fadeouts == 0){inits = abind(inits, inits_tmp, along = 1)}
  #else if(dim(inits)[1] == 999){inits = abind(inits, inits_tmp, along = 1)}
  else {inits = abind(inits, inits_tmp[-fadeouts_tmp,,], along = 1)}
  print(dim(inits))
  #if(dim(inits)[1] == 999){browser()}
}

# select ICs based on immunity
inits = inits[order(inits[,"R",t_unbounded + t_distanced]),,]
base_init = inits[500,,t_unbounded + t_distanced]

# saveRDS(list(full = inits, low = low_samp, mid = mid_samp, high = high_samp),output_initial_conditions_full)
saveRDS(list(full = inits),file.path(paste0(output_file,"/initial_sims_full.rda")))

# create initial conditions for prior immunity sensitivity
high_init = base_init
delta = 1000 - high_init["R"] 
high_init["R"] = high_init["R"] + delta

higher_init = base_init
delta = 2000 - high_init["R"] 
higher_init["S"] = higher_init["S"] - delta
higher_init["R"] = higher_init["R"] + delta


write.csv(t(base_init), file = file.path(paste0(output_file,"/initial_conditions_base.rda")), row.names = FALSE)
write.csv(t(high_init), file = file.path(paste0(output_file,"/initial_conditions_high.rda")), row.names = FALSE)
write.csv(t(higher_init), file = file.path(paste0(output_file,"/initial_conditions_higher.rda")), row.names = FALSE)