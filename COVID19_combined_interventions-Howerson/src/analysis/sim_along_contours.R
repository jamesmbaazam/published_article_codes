### Preamble
suppressWarnings({
  suppressMessages({
    require(tidyr)
    require(plyr)
    require(dplyr)
    require(pmultinom)
    require(abind)
    require(reshape)
    require(reshape2)
    require(metR)
    require(grid)
    require(wrswoR)
  })})



###### Read command line args ##########
args <- commandArgs(trailingOnly = TRUE)

input_model_sims <- args[1]
input_initial_conditions <- args[2]
input_params <- args[3]
input_model_file <- args[4]
input_model_tools <- args[5]
input_analysis_tools <- args[6]
input_nsims <- args[7]
output_contour_sims <- args[8]


# input_model_sims <- "output/sims/immediate/immediate_model_simulations.csv"
# input_initial_conditions <- "output/sims/initial_conditions/initial_conditions_base.csv"
# input_params <-"output/sims/parameters.rda"
# input_analysis_tools <- "src/analysis/analysis_tools.R"
# input_model_file <- "src/model/model_immediate_symp_isolation.R"
# input_model_tools <- "src/model/model_tools.R"
# input_analysis_tools <- "src/analysis/analysis_tools.R"
# input_nsims <- 5000
# output_contour_sims <- "output/sims/immediate/immediate_contours_simulations.csv"


source(input_model_file)
source(input_model_tools)
source(input_analysis_tools)
out = read.csv(input_model_sims)
inits = read.csv(input_initial_conditions)
params = readRDS(input_params)
list2env(params, globalenv())
#contour_lines = read.csv(input_contours)

sims = input_nsims

#### sim along contour ####

# calculate points along contour
out_sub = out %>% filter(objective == "tot.cases", tot.tests>0, Ts>0, sens == 1)


contour_lines = data.frame(tot.tests = integer(),
                           sens = double(),
                           level = integer(),
                           x = double(),
                           y = double())
for(i in unique(out_sub$tot.tests)){
    c = contour(out_sub %>% filter(tot.tests == i), "Ds", "d", "median", c(500))
    if(!is.null(c)){contour_lines = rbind(contour_lines,cbind(i,c))}
}
colnames(contour_lines) = c("tot.tests", "level", "Ds", "d")

# if multiple Ds for one d value (due to small wobbles in the contour from numerical approx), choose the smallest
contour_lines = contour_lines %>% group_by(tot.tests, level, d) %>% summarise(Ds = min(Ds))

tst.intervention = contour_lines %>% filter(tot.tests %in% c(100,500,2000,5000))
tst.intervention$contour_FLAG = 1 # is this on the contour - yes
tst.intervention2 = tst.intervention
tst.intervention2$d = tst.intervention2$d - 0.05
tst.intervention2$contour_FLAG = 0 # off contour (less distancing than necessary)
tst.intervention = rbind(tst.intervention, tst.intervention2)
tst.intervention = tst.intervention %>% filter(d >0)
tst.intervention$Da = tst.intervention$Ds
tst.intervention$Ts = 1/tst.intervention$Ds
tst.intervention$Ta = 1/tst.intervention$Da
tst.intervention$sim = 1:nrow(tst.intervention)
rm(tst.intervention2)

tst.asymp = data.frame(p = 0.4)
#tst.asymp = merge(tst.asymp, asymp[,c("p", "beta")])
#tst.asymp = merge(tst.asymp, start2)
tst.asymp$current.sim = 1:nrow(tst.asymp)

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

out = do.call(rbind, out)
out = do.call(rbind, out)
out$objective = rep(c("tot.cases", "peak.hosp", "days.above.thresh", "tot.hosp","death","sympiso","asympiso"), times = (nrow(out)/7))
colnames(out)[colnames(out) %in% c("quant.1","quant.2","quant.3")] = c("median", "lower","upper")

write.csv(out, output_contour_sims)

#saveRDS(out, "./FINAL CODE/output/sims/delayed_model_contours.rda")

# ggplot(data = tst.intervention)+
#   geom_point(aes(x = Ds, y = d, color = as.factor(tot.tests), shape = as.factor(contour_FLAG)))+
#   geom_line(aes(x = Ds, y = d, color = as.factor(tot.tests), group = paste(tot.tests, contour_FLAG)))
# 
