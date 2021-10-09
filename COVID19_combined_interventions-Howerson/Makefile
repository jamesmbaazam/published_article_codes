#######################
# Files and directories
# ---------------------
delayed_mod = "src/model/model_delayed_symp_isolation.R"
immed_mod = "src/model/model_immediate_symp_isolation.R"

model_tools = "src/model/model_tools.R"
calibration_tools = "src/model/calibration_tools.R"
analysis_tools = "src/analysis/analysis_tools.R"
viz_tools = "src/viz/viz_tools.R"

##################
# Necessary to do
# ----------------
parameters: 
	Rscript src/model/define_parameters.R \
	$(calibration_tools) \
	"output/sims/parameters.rda"

setup_sims: 
	Rscript src/analysis/setup_simulations.R \
		"output/sims/parameters.rda" \
		$(calibration_tools) \
		"output/sims/sims_setup.rda"
		

##################
# Model Simulation
# ----------------
# generate initial conditions for all model simulations
initial_conditions: parameters setup_sims
	Rscript src/analysis/generate_initial_conditions.R \
	$(immed_mod) \
	$(model_tools) \
	"output/sims/parameters.rda" \
	"output/sims/initial_conditions" \

# run simulations for each model
execute_immed_mod:
	Rscript src/analysis/simulate_model.R \
		$(immed_mod) \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		"output/sims/parameters.rda" \
		$(model_tools) \
		"output/sims/sims_setup.rda" \
		"main"\
		5000 \
		"output/sims/immediate/immediate_model_simulations.csv"
		
execute_delay_mod:
	Rscript src/analysis/simulate_model.R \
		$(delayed_mod) \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		"output/sims/parameters.rda" \
		$(model_tools) \
		"output/sims/sims_setup.rda" \
		"main"\
		5000 \
		"output/sims/delayed/delayed_model_simulations.csv"
		
# calculate interventions which yield equivalent outcomes
# simulate again along these contours
sim_contours_immed: 
	Rscript src/analysis/sim_along_contours.R \
		"output/sims/immediate/immediate_model_simulations.csv" \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		"output/sims/parameters.rda" \
		$(immed_mod) \
		$(model_tools) \
		$(analysis_tools) \
		5000 \
		"output/sims/immediate/immediate_contours_simulations.csv"
		
sim_contours_delay: 
	Rscript src/analysis/sim_along_contours.R \
		"output/sims/delayed/delayed_model_simulations.csv" \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		"output/sims/parameters.rda" \
		$(delayed_mod) \
		$(model_tools) \
		$(analysis_tools) \
		5000 \
		"output/sims/delayed/delayed_contours_simulations.csv"
		

#############
# Sensitivity
# -----------

execute_all_sens: execute_asymp_sens execute_asymp_rho_sens execute_IC_high_sens execute_IC_higher_sens

execute_all_asymp_sens: execute_asymp_sens execute_asymp_rho_sens

execute_all_IC_sens: execute_IC_high_sens execute_IC_higher_sens

execute_asymp_sens:
	Rscript src/analysis/simulate_model.R \
		$(immed_mod) \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		"output/sims/parameters.rda" \
		$(model_tools) \
		"output/sims/sims_setup.rda" \
		"sensitivity_asymp"\
		5000 \
		"output/sims/sensitivity/sensitivity_asymp_simulations.csv"

execute_asymp_rho_sens:
	Rscript src/analysis/simulate_model.R \
		$(immed_mod) \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		"output/sims/parameters.rda" \
		$(model_tools) \
		"output/sims/sims_setup.rda" \
		"sensitivity_asymp_rho"\
		5000 \
		"output/sims/sensitivity/sensitivity_asymp_rho_simulations.csv"
		

execute_IC_high_sens:
	Rscript src/analysis/simulate_model.R \
		$(immed_mod) \
		"output/sims/initial_conditions/initial_conditions_high.csv" \
		"output/sims/parameters.rda" \
		$(model_tools) \
		"output/sims/sims_setup.rda" \
		"sensitivity_IC"\
		5000 \
		"output/sims/sensitivity/sensitivity_IC_high_simulations.csv"
		
execute_IC_higher_sens:
	Rscript src/analysis/simulate_model.R \
		$(immed_mod) \
		"output/sims/initial_conditions/initial_conditions_higher.csv" \
		"output/sims/parameters.rda" \
		$(model_tools) \
		"output/sims/sims_setup.rda" \
		"sensitivity_IC"\
		5000 \
		"output/sims/sensitivity/sensitivity_IC_higher_simulations.csv"



###############
# Visualization
# -------------

all_figs: main_figs supp_figs

main_figs: Figure_2 Figure_3 Figure_4

supp_figs: Figure_S1 Figure_S2 Figure_S3 Figure_S4 Figure_S5 Figure_S6
		
# Figure 1 from Illustrator

Figure_2: 
	Rscript src/viz/heatmap_plots.R \
		"output/sims/immediate/immediate_model_simulations.csv" \
		$(analysis_tools) \
		$(viz_tools) \
		"TRUE" \
		"TRUE" \
		"FALSE" \
		"output/figures/Fig2.pdf"
		
Figure_3: 
	Rscript src/viz/contour_sims_plots.R \
		"output/sims/delayed/delayed_contours_simulations.csv" \
		"output/sims/immediate/immediate_contours_simulations.csv" \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		$(viz_tools) \
		"immed" \
		"TRUE" \
		"output/figures/Fig3.pdf"
		
Figure_4: 
	Rscript src/viz/change_outcomes_duncertainty.R \
		"output/sims/immediate/immediate_contours_simulations.csv" \
		$(viz_tools) \
		"output/figures/Fig4.pdf"

Figure_S1:
	Rscript src/viz/initial_conditions_plot.R \
		"output/sims/initial_conditions/initial_sims_full.rda" \
		"output/figures/FigS1.pdf"
		
Figure_S2:
	Rscript src/viz/heatmap_plots.R \
		"output/sims/immediate/immediate_model_simulations.csv" \
		$(analysis_tools) \
		$(viz_tools) \
		"FALSE" \
		"FALSE" \
		"TRUE" \
		"output/figures/FigS2.pdf"

Figure_S3: 
	Rscript src/viz/heatmap_plots.R \
		"output/sims/delayed/delayed_model_simulations.csv" \
		$(analysis_tools) \
		$(viz_tools) \
		"TRUE" \
		"FALSE" \
		"FALSE" \
		"output/figures/FigS3.pdf"

Figure_S4: 
	Rscript src/viz/sensitivity_plots.R \
		"output/sims/immediate/immediate_model_simulations.csv" \
		"output/sims/sensitivity" \
		$(viz_tools) \
		"output/figures/FigS4.pdf" 


Figure_S5: 
	Rscript src/viz/contour_sims_plots.R \
		"output/sims/delayed/delayed_contours_simulations.csv" \
		"output/sims/immediate/immediate_contours_simulations.csv" \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		$(viz_tools) \
		"diffiso" \
		"FALSE" \
		"output/figures/FigS5.pdf"
	
Figure_S6: 
	Rscript src/viz/contour_sims_plots.R \
		"output/sims/delayed/delayed_contours_simulations.csv" \
		"output/sims/immediate/immediate_contours_simulations.csv" \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		$(viz_tools) \
		"both" \
		"FALSE" \
		"output/figures/FigS6.pdf"
	
	
