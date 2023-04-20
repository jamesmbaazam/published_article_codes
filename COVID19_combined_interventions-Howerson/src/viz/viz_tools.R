# labels
labs_p = c("Asymptoamtic: 10%", "50%","90%")
names(labs_p) = c(0.1,0.5,0.9)
labs_p_notitle = c("10%", "50%","90%")
names(labs_p_notitle) = c(0.1,0.5,0.9)

labs_d = c("NPI Intensity: 0%","10%", "20%", "30%", "40%", "50%", "60%")
names(labs_d) =seq(0,0.6,0.1)
labs_d_notitle = c("0%","10%", "20%", "30%", "40%", "50%", "60%")
names(labs_d_notitle) =seq(0,0.6,0.1)

obj = c("tot.cases", "peak.hosp", "days.above.thresh", "tot.hosp", "death")
labs_obj = c("Infections","Peak Hospitalizations","Days Above Hospital Threshold", "Hospitalizations", "Deaths")
names(labs_obj) = obj
 
labs_tot_tests = c("Test administration: 1%", "5%", "20%", "50%")
names(labs_tot_tests) = c(100,500,2000,5000)

labs_Ds = c("Test delay: 1 hour", "1 day", "2 days", "5 days")
names(labs_Ds) = c(1/24,1,2,5)