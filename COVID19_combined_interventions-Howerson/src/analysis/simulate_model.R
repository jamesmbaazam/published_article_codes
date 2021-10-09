### Preamble
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
input_initial_conditions <- args[2]
input_params <- args[3]
input_model_tools <- args[4]
input_model_setup <- args[5]
input_sim_setup_type <- args[6]
input_nsims <- args[7]
output_sim_results_file <- args[8]


# input_model_file <-"src/model/model_immediate_symp_isolation.R"
# input_initial_conditions <-"output/sims/inital_conditions_mid.csv"
# input_params <- "output/sims/parameters.rda"
# input_model_tools <-"src/model/model_tools.R"
# input_model_setup <- "output/sims/sims_setup.rda"
# input_sim_setup_type <- "main"
# output_sim_results_file <- "output/sims/immediate_model_simulations.csv"

# input_model_file <-"src/model/model_immediate_symp_isolation.R"
# input_initial_conditions <- "output/sims/inital_conditions_mid.csv" 
# input_params <- "output/sims/parameters.rda"
# input_model_tools <-"src/model/model_tools.R"
# input_model_setup <- "output/sims/sims_setup.rda"
# input_sim_setup_type <- "sensitivity_asymp"
# output_sim_results_file <- "output/sims/sensitivity_asymp_simulations.csv"

##### simulate model #####

source(input_model_file)
source(input_model_tools)
inits = read.csv(input_initial_conditions)
params = readRDS(input_params)
list2env(params, globalenv())
setup = readRDS(input_model_setup)

sims = input_nsims

tst.intervention = setup[[input_sim_setup_type]][["int"]]
tst.asymp = setup[[input_sim_setup_type]][["asymp"]]

print(paste("starting", sims, "simulations"))
print(paste("for", nrow(tst.intervention), "interventions"))

start_time = Sys.time()
out = apply(tst.asymp,1,fixedtest_apply, # note: inits = tst.asymp
            change.vars = tst.intervention, 
            params = params, 
            times.first = t_distanced-1, 
            times.second = 30, 
            IC = fixedtest_setup_IC(inits,sims), 
            sims = sims) 
end_time = Sys.time()
end_time-start_time

out = do.call(rbind,out)
out = do.call(rbind,out)
out$objective = rep(c("tot.cases", "peak.hosp", "days.above.thresh", "tot.hosp","death", "sympiso", "asympiso"), times = (nrow(out)/7))
colnames(out)[colnames(out) %in% c("quant.1","quant.2","quant.3")] = c("median", "lower","upper")


#saveRDS(out, output_sim_results_file)

write.csv(out, output_sim_results_file, row.names = FALSE)
