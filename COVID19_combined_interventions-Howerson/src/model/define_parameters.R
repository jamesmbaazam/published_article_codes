###### Read command line args ##########
args <- commandArgs(trailingOnly = TRUE)

input_calibration_tools <- args[1]
output_params_file <- args[2]

# input_calibration_tools <- "src/model/calibration_tools.R"
# output_params_file <- "output/sims/parameters.rda"

###### Define parameters for analysis ####

# load calibration functions
source(input_calibration_tools)


#sims = 5000 # number of simulations to perform for each intervention

t_unbounded = 10 # number of days epidemic grows unbounded
t_distanced = 90 # number of days epidemic grows constrained
d_distanced = 0.55 # d during distanced period

# define parameters to vary in analysis
test_d = seq(0,0.6, by = 0.05)
test_Ds = c(1/24,seq(0.5,8,by = 0.5))
test_tottests = c(0,100,250,500,1000,1500,2000,2500,3000,3500,4000,4500,5000)
test_sens = c(1,0.9)




#### set up parameters ####
## biological parameters (for pre-intervention epidemic)
rho = 1 # proportion of infectiousness of asymptoamtic infections (betaA = rho*betaI), will consider 0.5 in sens - from Li et al. Science 2020. 
# model parameters from Davies et al. 
theta = 1/3 # exposed to infected
gamma_I1I2 = 1/2.1 # presymptoamtic to reported symptoms
gamma_AR = 1/5 # asymptoamtic recovery rate
gamma_I2R = 1/2.9 # mild symptoamtic recovery rate
gamma_I4H = 1/2.1 # time from severe symptom reporting to hospitalization, CHOSEN FOR EQUIVALENT TIMING ACROSS TRACKS
gamma_HTR = 1/12 # hospitalization to recovery
p = 0.4 # proportion who never report symptoms (sub-clinical) - see reviews by Oran and Topol and Buitrago-Garcia et al. 
q = 0.05 # proportion of symptomatic cases that are hospitalized
gamma_HD = 1/7.5 # death rate given hospitalization,  
alpha = 1/7.5 # another death rate param... delete 
c = 0.001 # background test seeking

H.thresh = 50 # hospitalization threshold

## baseline intervention params
# testing
tot.test = 0 # no tests at beginning of outbreak
sens = 1 # true positive rate (1-sens: false negative)
spec = 1 # true negative rate (1-spec: false positive)
Ta = 0 # asymptomatic isolation rate
Ts = 0 # symptomatic isolation rate
# distancing
d = 0 # no distancing at beginning of outbreak


params = c(rho = rho, theta = theta, d = d, p = p,q = q,Ta = Ta,Ts = Ts, gamma_HD = gamma_HD,alpha = alpha,
           gamma_I1I2 = gamma_I1I2, gamma_I2R = gamma_I2R, gamma_AR = gamma_AR, gamma_I2H = gamma_I4H, gamma_I4H = gamma_I4H, gamma_HTR = gamma_HTR,
           tot.tests = tot.test, sens = sens, spec = spec, c = c) ## FOR NOW - REMOVE BETA AND CALC

# R0 between 2 and 3, fit beta s.t. R0 = 2.5

#### calc beta - best guess p ####
# configure model for NGM ####
mod = mod_config()
# find beta s.t. R0 = 2.5
rslts = data.frame(beta = seq(0.4,0.6, by = 0.001))
rslts = apply(rslts, 1, apply.nextGen, in.change.params.names = c("beta"),
              in.istates = mod$istates, in.flist = mod$flist, in.vlist = mod$vlist, 
              in.params = params, in.dfe = mod$df, in.sim = 1)
rslts = do.call(rbind, rslts)
est.beta = data.frame(R0 = c(1,2,2.5,3), beta = NA)
for(i in est.beta$R0){
  est.beta$beta[which(est.beta$R0 == i)] = rslts$beta[which(abs(rslts$R0-i) == min(abs(rslts$R0-i)))]
}
# beta s.t. R0 = 2.5
#est.beta$beta[est.beta$R0 == 2.5]
params["beta"] = est.beta$beta[est.beta$R0 == 2.5]


saveRDS(list(params = params,t_unbounded = t_unbounded, t_distanced = t_distanced, d_distanced = d_distanced,
             test_d = test_d, test_Ds = test_Ds, test_sens = test_sens, test_tottests = test_tottests, H.thresh = H.thresh), output_params_file) #sims = sims, 