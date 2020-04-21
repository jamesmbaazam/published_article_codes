####This file contains the code to run the 37 models with and wihout interventions in the 
####paper titled "Essential information: Uncertainty and optimal control of Ebola outbreaks"
#### Before running the code, you would need to first import the files named Functions.R and Parameters.R
#Import project-specific functions and model parameters
source('Functions.R')
source('Parameters.R')

nrun = 100;			# number of simulation runs
times=seq(0,2000,by=0.25)	# time steps
save.var=c("times","Sc","Ec","Ic1","Ih2","CC","CD","CP") 	# variables to save 
intensity <- seq(0.1,1,by = 0.1);		# The degree to which the intervention reduces (or increase for hospitalization) the relevant rate in the disease model

## Run each of the 37 models without interventions
runModel(pars.agu15,'agu15',nrun);
runModel(pars.cam14,'cam14',nrun);
runModel(pars.eli14_full_all,'eli14_full_all',nrun);
runModel(pars.gom14_SEIHFR,'gom14_SEIHFR',nrun);
runModel(pars.leg07_con95,'leg07_con95',nrun);
runModel(pars.leg07_uga00,'leg07_uga00',nrun);
runModel(pars.riv14_lib,'riv14_lib',nrun);
runModel(pars.riv14_sie,'riv14_sie',nrun);
runModel(pars.cho15,'cho15',nrun);
runModel(pars.fas14,'fas14',nrun);
runModel(pars.kha14_lib,'kha14_lib',nrun);
runModel(pars.kha14_sie,'kha14_sie',nrun);
runModel(pars.kuc15,'kuc15',nrun);
runModel(pars.eli14_part_all,'eli14_part_all',nrun);
runModel(pars.eli14_part_gui,'eli14_part_gui',nrun);
runModel(pars.eli14_part_lib,'eli14_part_lib',nrun);
runModel(pars.eli14_part_sie,'eli14_part_sie',nrun);
runModel(pars.wei15_gui,'wei15_gui',nrun);
runModel(pars.wei15_lib,'wei15_lib',nrun);
runModel(pars.wei15_sie,'wei15_sie',nrun);
runModel(pars.alt14_gui,'alt14_gui',nrun);
runModel(pars.alt14_lib,'alt14_lib',nrun);
runModel(pars.alt14_sie,'alt14_sie',nrun);
runModel(pars.bas14_gui,'bas14_gui',nrun);
runModel(pars.bas14_sie,'bas14_sie',nrun);
runModel(pars.cho04_con95,'cho04_con95',nrun);
runModel(pars.cho04_uga00,'cho04_uga00',nrun);
runModel(pars.gom14_SEIR,'gom14_SEIR',nrun);
runModel(pars.lek06_vag,'lek06_vag',nrun);
runModel(pars.lek06_inf,'lek06_inf',nrun);
runModel(pars.mel14_lib,'mel14_lib',nrun);
runModel(pars.mel14_sie,'mel14_sie',nrun);
runModel(pars.fer04_con95,'fer04_con95',nrun);
runModel(pars.fer04_uga00,'fer04_uga00',nrun);
runModel(pars.sha14_gui,'sha14_gui',nrun);
runModel(pars.sha14_lib,'sha14_lib',nrun);
runModel(pars.sha14_sie,'sha14_sie',nrun);

##Run each of the 37 models under the intervention of reducing funeral transmission
#full SEIHFR models
runModel.funtra.SEIHFR(pars.agu15,'agu15',nrun,intensity);
runModel.funtra.SEIHFR(pars.cam14,'cam14',nrun,intensity);
runModel.funtra.SEIHFR(pars.eli14_full_all,'eli14_full_all',nrun,intensity);
runModel.funtra.SEIHFR(pars.gom14_SEIHFR,'gom14_SEIHFR',nrun,intensity);
runModel.funtra.SEIHFR(pars.leg07_con95,'leg07_con95',nrun,intensity);
runModel.funtra.SEIHFR(pars.leg07_uga00,'leg07_uga00',nrun,intensity);
runModel.funtra.SEIHFR(pars.riv14_lib,'riv14_lib',nrun,intensity);
runModel.funtra.SEIHFR(pars.riv14_sie,'riv14_sie',nrun,intensity);
#SEIHR models
runModel.funtra.SEIHR(pars.cho15,'cho15',nrun,intensity);
runModel.funtra.SEIHR(pars.fas14,'fas14',nrun,intensity);
runModel.funtra.SEIHR(pars.kha14_lib,'kha14_lib',nrun,intensity);
runModel.funtra.SEIHR(pars.kha14_sie,'kha14_sie',nrun,intensity);
runModel.funtra.SEIHR(pars.kuc15,'kuc15',nrun,intensity);
#SEIFR models  
runModel.funtra.SEIFR(pars.eli14_part_all,'eli14_part_all',nrun,intensity);
runModel.funtra.SEIFR(pars.eli14_part_gui,'eli14_part_gui',nrun,intensity);
runModel.funtra.SEIFR(pars.eli14_part_lib,'eli14_part_lib',nrun,intensity);
runModel.funtra.SEIFR(pars.eli14_part_sie,'eli14_part_sie',nrun,intensity);
runModel.funtra.SEIFR(pars.wei15_gui,'wei15_gui',nrun,intensity);
runModel.funtra.SEIFR(pars.wei15_lib,'wei15_lib',nrun,intensity);
runModel.funtra.SEIFR(pars.wei15_sie,'wei15_sie',nrun,intensity);
#SEIR models
runModel.funtra.SEIR(pars.alt14_gui,'alt14_gui',nrun,intensity);
runModel.funtra.SEIR(pars.alt14_lib,'alt14_lib',nrun,intensity);
runModel.funtra.SEIR(pars.alt14_sie,'alt14_sie',nrun,intensity);
runModel.funtra.SEIR(pars.bas14_gui,'bas14_gui',nrun,intensity);
runModel.funtra.SEIR(pars.bas14_sie,'bas14_sie',nrun,intensity);
runModel.funtra.SEIR(pars.cho04_con95,'cho04_con95',nrun,intensity);
runModel.funtra.SEIR(pars.cho04_uga00,'cho04_uga00',nrun,intensity);
runModel.funtra.SEIR(pars.gom14_SEIR,'gom14_SEIR',nrun,intensity);
runModel.funtra.SEIR(pars.lek06_vag,'lek06_vag',nrun,intensity);
runModel.funtra.SEIR(pars.lek06_inf,'lek06_inf',nrun,intensity);
runModel.funtra.SEIR(pars.mel14_lib,'mel14_lib',nrun,intensity);
runModel.funtra.SEIR(pars.mel14_sie,'mel14_sie',nrun,intensity);
runModel.funtra.SEIR(pars.fer04_con95,'fer04_con95',nrun,intensity);
runModel.funtra.SEIR(pars.fer04_uga00,'fer04_uga00',nrun,intensity);
runModel.funtra.SEIR(pars.sha14_gui,'sha14_gui',nrun,intensity);
runModel.funtra.SEIR(pars.sha14_lib,'sha14_lib',nrun,intensity);
runModel.funtra.SEIR(pars.sha14_sie,'sha14_sie',nrun,intensity);

##Run each of the 37 models under the intervention of reducing community transmission
#full SEIHFR models
runModel.comtra.SEIHFR(pars.agu15,'agu15',nrun,intensity);
runModel.comtra.SEIHFR(pars.cam14,'cam14',nrun,intensity);
runModel.comtra.SEIHFR(pars.eli14_full_all,'eli14_full_all',nrun,intensity);
runModel.comtra.SEIHFR(pars.gom14_SEIHFR,'gom14_SEIHFR',nrun,intensity);
runModel.comtra.SEIHFR(pars.leg07_con95,'leg07_con95',nrun,intensity);
runModel.comtra.SEIHFR(pars.leg07_uga00,'leg07_uga00',nrun,intensity);
runModel.comtra.SEIHFR(pars.riv14_lib,'riv14_lib',nrun,intensity);
runModel.comtra.SEIHFR(pars.riv14_sie,'riv14_sie',nrun,intensity);
#SEIHR models
runModel.comtra.SEIHR(pars.cho15,'cho15',nrun,intensity);
runModel.comtra.SEIHR(pars.fas14,'fas14',nrun,intensity);
runModel.comtra.SEIHR(pars.kha14_lib,'kha14_lib',nrun,intensity);
runModel.comtra.SEIHR(pars.kha14_sie,'kha14_sie',nrun,intensity);
runModel.comtra.SEIHR(pars.kuc15,'kuc15',nrun,intensity);
#SEIFR models  
runModel.comtra.SEIFR(pars.eli14_part_all,'eli14_part_all',nrun,intensity);
runModel.comtra.SEIFR(pars.eli14_part_gui,'eli14_part_gui',nrun,intensity);
runModel.comtra.SEIFR(pars.eli14_part_lib,'eli14_part_lib',nrun,intensity);
runModel.comtra.SEIFR(pars.eli14_part_sie,'eli14_part_sie',nrun,intensity);
runModel.comtra.SEIFR(pars.wei15_gui,'wei15_gui',nrun,intensity);
runModel.comtra.SEIFR(pars.wei15_lib,'wei15_lib',nrun,intensity);
runModel.comtra.SEIFR(pars.wei15_sie,'wei15_sie',nrun,intensity);
#SEIR models
runModel.comtra.SEIR(pars.alt14_gui,'alt14_gui',nrun,intensity);
runModel.comtra.SEIR(pars.alt14_lib,'alt14_lib',nrun,intensity);
runModel.comtra.SEIR(pars.alt14_sie,'alt14_sie',nrun,intensity);
runModel.comtra.SEIR(pars.bas14_gui,'bas14_gui',nrun,intensity);
runModel.comtra.SEIR(pars.bas14_sie,'bas14_sie',nrun,intensity);
runModel.comtra.SEIR(pars.cho04_con95,'cho04_con95',nrun,intensity);
runModel.comtra.SEIR(pars.cho04_uga00,'cho04_uga00',nrun,intensity);
runModel.comtra.SEIR(pars.gom14_SEIR,'gom14_SEIR',nrun,intensity);
runModel.comtra.SEIR(pars.lek06_vag,'lek06_vag',nrun,intensity);
runModel.comtra.SEIR(pars.lek06_inf,'lek06_inf',nrun,intensity);
runModel.comtra.SEIR(pars.mel14_lib,'mel14_lib',nrun,intensity);
runModel.comtra.SEIR(pars.mel14_sie,'mel14_sie',nrun,intensity);
runModel.comtra.SEIR(pars.fer04_con95,'fer04_con95',nrun,intensity);
runModel.comtra.SEIR(pars.fer04_uga00,'fer04_uga00',nrun,intensity);
runModel.comtra.SEIR(pars.sha14_gui,'sha14_gui',nrun,intensity);
runModel.comtra.SEIR(pars.sha14_lib,'sha14_lib',nrun,intensity);
runModel.comtra.SEIR(pars.sha14_sie,'sha14_sie',nrun,intensity);

##Run each of the 37 models under the intervention of reducing mortality ratio
#full SEIHFR models
runModel.mor.SEIHFR(pars.agu15,'agu15',nrun,intensity);
runModel.mor.SEIHFR(pars.cam14,'cam14',nrun,intensity);
runModel.mor.SEIHFR(pars.eli14_full_all,'eli14_full_all',nrun,intensity);
runModel.mor.SEIHFR(pars.gom14_SEIHFR,'gom14_SEIHFR',nrun,intensity);
runModel.mor.SEIHFR(pars.leg07_con95,'leg07_con95',nrun,intensity);
runModel.mor.SEIHFR(pars.leg07_uga00,'leg07_uga00',nrun,intensity);
runModel.mor.SEIHFR(pars.riv14_lib,'riv14_lib',nrun,intensity);
runModel.mor.SEIHFR(pars.riv14_sie,'riv14_sie',nrun,intensity);
#SEIHR models
runModel.mor.SEIHR(pars.cho15,'cho15',nrun,intensity);
runModel.mor.SEIHR(pars.fas14,'fas14',nrun,intensity);
runModel.mor.SEIHR(pars.kha14_lib,'kha14_lib',nrun,intensity);
runModel.mor.SEIHR(pars.kha14_sie,'kha14_sie',nrun,intensity);
runModel.mor.SEIHR(pars.kuc15,'kuc15',nrun,intensity);
#SEIFR models  
runModel.mor.SEIFR(pars.eli14_part_all,'eli14_part_all',nrun,intensity);
runModel.mor.SEIFR(pars.eli14_part_gui,'eli14_part_gui',nrun,intensity);
runModel.mor.SEIFR(pars.eli14_part_lib,'eli14_part_lib',nrun,intensity);
runModel.mor.SEIFR(pars.eli14_part_sie,'eli14_part_sie',nrun,intensity);
runModel.mor.SEIFR(pars.wei15_gui,'wei15_gui',nrun,intensity);
runModel.mor.SEIFR(pars.wei15_lib,'wei15_lib',nrun,intensity);
runModel.mor.SEIFR(pars.wei15_sie,'wei15_sie',nrun,intensity);
#SEIR models
runModel.mor.SEIR(pars.alt14_gui,'alt14_gui',nrun,intensity);
runModel.mor.SEIR(pars.alt14_lib,'alt14_lib',nrun,intensity);
runModel.mor.SEIR(pars.alt14_sie,'alt14_sie',nrun,intensity);
runModel.mor.SEIR(pars.bas14_gui,'bas14_gui',nrun,intensity);
runModel.mor.SEIR(pars.bas14_sie,'bas14_sie',nrun,intensity);
runModel.mor.SEIR(pars.cho04_con95,'cho04_con95',nrun,intensity);
runModel.mor.SEIR(pars.cho04_uga00,'cho04_uga00',nrun,intensity);
runModel.mor.SEIR(pars.gom14_SEIR,'gom14_SEIR',nrun,intensity);
runModel.mor.SEIR(pars.lek06_vag,'lek06_vag',nrun,intensity);
runModel.mor.SEIR(pars.lek06_inf,'lek06_inf',nrun,intensity);
runModel.mor.SEIR(pars.mel14_lib,'mel14_lib',nrun,intensity);
runModel.mor.SEIR(pars.mel14_sie,'mel14_sie',nrun,intensity);
runModel.mor.SEIR(pars.fer04_con95,'fer04_con95',nrun,intensity);
runModel.mor.SEIR(pars.fer04_uga00,'fer04_uga00',nrun,intensity);
runModel.mor.SEIR(pars.sha14_gui,'sha14_gui',nrun,intensity);
runModel.mor.SEIR(pars.sha14_lib,'sha14_lib',nrun,intensity);
runModel.mor.SEIR(pars.sha14_sie,'sha14_sie',nrun,intensity);

##Run each of the 37 models under the intervention of reducing hospital transmission
#full SEIHFR models
runModel.hostra.SEIHFR(pars.agu15,'agu15',nrun,intensity);
runModel.hostra.SEIHFR(pars.cam14,'cam14',nrun,intensity);
runModel.hostra.SEIHFR(pars.eli14_full_all,'eli14_full_all',nrun,intensity);
runModel.hostra.SEIHFR(pars.gom14_SEIHFR,'gom14_SEIHFR',nrun,intensity);
runModel.hostra.SEIHFR(pars.leg07_con95,'leg07_con95',nrun,intensity);
runModel.hostra.SEIHFR(pars.leg07_uga00,'leg07_uga00',nrun,intensity);
runModel.hostra.SEIHFR(pars.riv14_lib,'riv14_lib',nrun,intensity);
runModel.hostra.SEIHFR(pars.riv14_sie,'riv14_sie',nrun,intensity);
#SEIHR models
runModel.hostra.SEIHR(pars.cho15,'cho15',nrun,intensity);
runModel.hostra.SEIHR(pars.fas14,'fas14',nrun,intensity);
runModel.hostra.SEIHR(pars.kha14_lib,'kha14_lib',nrun,intensity);
runModel.hostra.SEIHR(pars.kha14_sie,'kha14_sie',nrun,intensity);
runModel.hostra.SEIHR(pars.kuc15,'kuc15',nrun,intensity);
#SEIFR models  
runModel.hostra.SEIFR(pars.eli14_part_all,'eli14_part_all',nrun,intensity);
runModel.hostra.SEIFR(pars.eli14_part_gui,'eli14_part_gui',nrun,intensity);
runModel.hostra.SEIFR(pars.eli14_part_lib,'eli14_part_lib',nrun,intensity);
runModel.hostra.SEIFR(pars.eli14_part_sie,'eli14_part_sie',nrun,intensity);
runModel.hostra.SEIFR(pars.wei15_gui,'wei15_gui',nrun,intensity);
runModel.hostra.SEIFR(pars.wei15_lib,'wei15_lib',nrun,intensity);
runModel.hostra.SEIFR(pars.wei15_sie,'wei15_sie',nrun,intensity);
#SEIR models
runModel.hostra.SEIR(pars.alt14_gui,'alt14_gui',nrun,intensity);
runModel.hostra.SEIR(pars.alt14_lib,'alt14_lib',nrun,intensity);
runModel.hostra.SEIR(pars.alt14_sie,'alt14_sie',nrun,intensity);
runModel.hostra.SEIR(pars.bas14_gui,'bas14_gui',nrun,intensity);
runModel.hostra.SEIR(pars.bas14_sie,'bas14_sie',nrun,intensity);
runModel.hostra.SEIR(pars.cho04_con95,'cho04_con95',nrun,intensity);
runModel.hostra.SEIR(pars.cho04_uga00,'cho04_uga00',nrun,intensity);
runModel.hostra.SEIR(pars.gom14_SEIR,'gom14_SEIR',nrun,intensity);
runModel.hostra.SEIR(pars.lek06_vag,'lek06_vag',nrun,intensity);
runModel.hostra.SEIR(pars.lek06_inf,'lek06_inf',nrun,intensity);
runModel.hostra.SEIR(pars.mel14_lib,'mel14_lib',nrun,intensity);
runModel.hostra.SEIR(pars.mel14_sie,'mel14_sie',nrun,intensity);
runModel.hostra.SEIR(pars.fer04_con95,'fer04_con95',nrun,intensity);
runModel.hostra.SEIR(pars.fer04_uga00,'fer04_uga00',nrun,intensity);
runModel.hostra.SEIR(pars.sha14_gui,'sha14_gui',nrun,intensity);
runModel.hostra.SEIR(pars.sha14_lib,'sha14_lib',nrun,intensity);
runModel.hostra.SEIR(pars.sha14_sie,'sha14_sie',nrun,intensity);

##Run each of the 37 models under the intervention of increasing hospitalization
#full SEIHFR models
runModel.hospitalization.SEIHFR(pars.agu15,'agu15',nrun,intensity);
runModel.hospitalization.SEIHFR(pars.cam14,'cam14',nrun,intensity);
runModel.hospitalization.SEIHFR(pars.eli14_full_all,'eli14_full_all',nrun,intensity);
runModel.hospitalization.SEIHFR(pars.gom14_SEIHFR,'gom14_SEIHFR',nrun,intensity);
runModel.hospitalization.SEIHFR(pars.leg07_con95,'leg07_con95',nrun,intensity);
runModel.hospitalization.SEIHFR(pars.leg07_uga00,'leg07_uga00',nrun,intensity);
runModel.hospitalization.SEIHFR(pars.riv14_lib,'riv14_lib',nrun,intensity);
runModel.hospitalization.SEIHFR(pars.riv14_sie,'riv14_sie',nrun,intensity);
#SEIHR models
runModel.hospitalization.SEIHR(pars.cho15,'cho15',nrun,intensity);
runModel.hospitalization.SEIHR(pars.fas14,'fas14',nrun,intensity);
runModel.hospitalization.SEIHR(pars.kha14_lib,'kha14_lib',nrun,intensity);
runModel.hospitalization.SEIHR(pars.kha14_sie,'kha14_sie',nrun,intensity);
runModel.hospitalization.SEIHR(pars.kuc15,'kuc15',nrun,intensity);
#SEIFR models  
runModel.hospitalization.SEIFR(pars.eli14_part_all,'eli14_part_all',nrun,intensity);
runModel.hospitalization.SEIFR(pars.eli14_part_gui,'eli14_part_gui',nrun,intensity);
runModel.hospitalization.SEIFR(pars.eli14_part_lib,'eli14_part_lib',nrun,intensity);
runModel.hospitalization.SEIFR(pars.eli14_part_sie,'eli14_part_sie',nrun,intensity);
runModel.hospitalization.SEIFR(pars.wei15_gui,'wei15_gui',nrun,intensity);
runModel.hospitalization.SEIFR(pars.wei15_lib,'wei15_lib',nrun,intensity);
runModel.hospitalization.SEIFR(pars.wei15_sie,'wei15_sie',nrun,intensity);
#SEIR models
runModel.hospitalization.SEIR(pars.alt14_gui,'alt14_gui',nrun,intensity);
runModel.hospitalization.SEIR(pars.alt14_lib,'alt14_lib',nrun,intensity);
runModel.hospitalization.SEIR(pars.alt14_sie,'alt14_sie',nrun,intensity);
runModel.hospitalization.SEIR(pars.bas14_gui,'bas14_gui',nrun,intensity);
runModel.hospitalization.SEIR(pars.bas14_sie,'bas14_sie',nrun,intensity);
runModel.hospitalization.SEIR(pars.cho04_con95,'cho04_con95',nrun,intensity);
runModel.hospitalization.SEIR(pars.cho04_uga00,'cho04_uga00',nrun,intensity);
runModel.hospitalization.SEIR(pars.gom14_SEIR,'gom14_SEIR',nrun,intensity);
runModel.hospitalization.SEIR(pars.lek06_vag,'lek06_vag',nrun,intensity);
runModel.hospitalization.SEIR(pars.lek06_inf,'lek06_inf',nrun,intensity);
runModel.hospitalization.SEIR(pars.mel14_lib,'mel14_lib',nrun,intensity);
runModel.hospitalization.SEIR(pars.mel14_sie,'mel14_sie',nrun,intensity);
runModel.hospitalization.SEIR(pars.fer04_con95,'fer04_con95',nrun,intensity);
runModel.hospitalization.SEIR(pars.fer04_uga00,'fer04_uga00',nrun,intensity);
runModel.hospitalization.SEIR(pars.sha14_gui,'sha14_gui',nrun,intensity);
runModel.hospitalization.SEIR(pars.sha14_lib,'sha14_lib',nrun,intensity);
runModel.hospitalization.SEIR(pars.sha14_sie,'sha14_sie',nrun,intensity);







