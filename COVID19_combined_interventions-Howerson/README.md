# Synergistic  interventions to control COVID-19: mass testing and isolation mitigates reliance on distancing
This repository contains code to reproduce results in Howerton et al. *Submitted*. All code is licensed under the CC-BY-NC Creative Commons attribution-noncommercial license (http://creativecommons.org/licenses/by-nc/3.0/).

## Code structure
    code
    ├── LICENSE
    ├── Makefile
    ├── src
    │   ├── analysis
    │   │   ├── analysis_tools.R
    │   │   ├── generate_initial_conditions.R
    │   │   ├── setup_simulations.R
    │   │   ├── sim_along_contours.R
    |   |   └── simulate_model.R
    │   ├── model
    │   │   ├── calibration_tools.R
    │   │   ├── define_parameters.R
    │   │   ├── model_delayed_symp_isolation.R
    │   │   ├── model_immediate_symp_isolation.R
    |   |   └── model_tools.R  
    │   ├── viz
    │   │   ├── change_outcomes_duncertainty.R
    │   │   ├── contour_sims_plots.R
    │   │   ├── correlation_plot.R
    │   │   ├── heatmap_plots.R
    │   │   ├── initial_conditions_plot.R
    │   │   ├── sensitivity_plots.R
    |   └── └── viz_tools.R  
    ├── output
    │   ├── figures
    │   ├── sims
    │   │   ├── delayed
    │   │   ├── immediate
    │   │   ├── initial_conditions
    └── └── └── sensitivity

## Requirements

This analysis was run in R >3.6 and uses several external packages.  All required R packages are listed in [`src/analysis/simulate_model.R`](src/analysis/simulate_model.R).  The analysis can be controlled using GNU Make (see the [Makefile](Makefile) for details) or each subanalysis run using the commands listed below.  


## Analysis 
Simulation data can be generated using the following steps: 

1. **Setup:** Configure model and generate initial conditions
2. **Simulate interventions:** Run model simulations under varying intervention scenarios
3. **Simulate along contours:** Identify interventions with equivalent outcomes, and resimulate these

Data used for manuscript is available in `output\sims`.  Scripts can also be executed using `make` command via `Makefile`.

### 1. Configure model and generate initial conditions

*Step 1:* create .rda object with model parameters (including biological and intervention)
```bash
Rscript src/model/define_parameters.R \
	"src/model/calibration_tools.R" \
	"output/sims/parameters.rda"
```
*Step 2:* create second .rda object that defines set of intervention combinations to simulate.
```bash
Rscript src/analysis/setup_simulations.R \
	"output/sims/parameters.rda" \
	"src/model/calibration_tools.R" \
	"output/sims/sims_setup.rda"
```
*Step 3:* generate initial conditions for interventions
```bash
Rscript src/analysis/generate_initial_conditions.R \
	"src/model/model_immediate_symp_isolation.R" \
	"src/model/model_tools.R" \
	"output/sims/parameters.rda" \
	"output/sims/initial_conditions" \
```

### 2. Simulate interventions
*Step 1:* model assuming symptomatic infections quarantine immediately upon test administration, whereas asymptomatic infections isolate only upon positive test
```bash
Rscript src/analysis/simulate_model.R \
	"src/model/model_immediate_symp_isolation.R" \
	"output/sims/initial_conditions/initial_conditions_base.csv" \
	"output/sims/parameters.rda" \
	"src/model/model_tools.R" \
	"output/sims/sims_setup.rda" \
	"main"\
	5000 \
	"output/sims/immediate/immediate_model_simulations.csv"
```
*Step 2:* model assuming all infections wait for a positive test result to isolate
```bash
Rscript src/analysis/simulate_model.R \
	"src/model/model_delayed_symp_isolation.R" \
	"output/sims/initial_conditions/initial_conditions_base.csv" \
	"output/sims/parameters.rda" \
	"src/model/model_tools.R" \
	"output/sims/sims_setup.rda" \
	"main"\
	5000 \
	"output/sims/delayed/delayed_model_simulations.csv"
```

### 3. Simulations along contours 
*Step 1:* execute for immediate model
```bash
Rscript src/analysis/sim_along_contours.R \
	"output/sims/immediate/immediate_model_simulations.csv" \
	"output/sims/initial_conditions/initial_conditions_base.csv" \
	"output/sims/parameters.rda" \
	"src/model/model_immediate_symp_isolation.R" \
	"src/model/model_tools.R" \
	"src/analysis/analysis_tools.R" \
	5000 \
	"output/sims/immediate/immediate_contours_simulations.csv"
```
*Step 2:* execute for delayed model
```bash
Rscript src/analysis/sim_along_contours.R \
	"output/sims/delayed/delayed_model_simulations.csv" \
	"output/sims/initial_conditions/initial_conditions_base.csv" \
	"output/sims/parameters.rda" \
	"src/model/model_delayed_symp_isolation.R" \
	"src/model/model_tools.R" \
	"src/analysis/analysis_tools.R" \
	5000 \
	"output/sims/delayed/delayed_contours_simulations.csv"
```

## Figures

### Figure 1
Created in Adobe illustrator, no code available. 

![./output/figures/Fig1.pdf](./output/figures/Fig1.pdf)

### Figure 2
```bash

Rscript src/viz/heatmap_plots.R \
		"output/sims/immediate/immediate_model_simulations.csv" \
		"src/analysis/analysis_tools.R" \
		"src/viz/viz_tools.R" \
		"TRUE" \
		"TRUE" \
		"FALSE" \
		"output/figures/Fig2.pdf"
```

![./output/figures/Fig2.pdf](./output/figures/Fig2.pdf)


### Figure 3
Note availability labels added in Adobe Illustrator.
```bash
Rscript src/viz/contour_sims_plots.R \
		"output/sims/delayed/delayed_contours_simulations.csv" \
		"output/sims/immediate/immediate_contours_simulations.csv" \
		"output/sims/initial_conditions/initial_conditions_base.csv" \
		"src/viz/viz_tools.R" \
		"immed" \
		"TRUE" \
		"output/figures/Fig3.pdf"
```

![./output/figures/Fig3.pdf](./output/figures/Fig3.pdf)


### Figure 4
```bash
Rscript src/viz/change_outcomes_duncertainty.R \
	"output/sims/immediate/immediate_contours_simulations.csv" \
	"src/viz/viz_tools.R" \
	"output/figures/Fig4.pdf"
```

![./output/figures/Fig4.pdf](./output/figures/Fig4.pdf)

### Supplementary Figure 1
```bash
Rscript src/viz/initial_conditions_plot.R \
		"output/sims/initial_conditions/initial_sims_full.rda" \
		"output/figures/FigS1.pdf"
```

![./output/figures/FigS1.pdf](./output/figures/FigS1.pdf)


### Supplementary Figure 2
```bash
Rscript src/viz/heatmap_plots.R \
	"output/sims/immediate/immediate_model_simulations.csv" \
	$(analysis_tools) \
	$(viz_tools) \
	"FALSE" \
	"FALSE" \
	"TRUE" \
	"output/figures/FigS2.pdf"
```

![./output/figures/FigS2.pdf](./output/figures/FigS2.pdf)


### Supplementary Figure 3
```bash
Rscript src/viz/heatmap_plots.R \
	"output/sims/delayed/delayed_model_simulations.csv" \
	$(analysis_tools) \
	$(viz_tools) \
	"TRUE" \
	"FALSE" \
	"FALSE" \
	"output/figures/FigS3.pdf"
```

![./output/figures/FigS3.pdf](./output/figures/FigS3.pdf)


### Supplementary Figure 4
```bash
Rscript src/viz/sensitivity_plots.R \
	"output/sims/immediate/immediate_model_simulations.csv" \
	"output/sims/sensitivity" \
	$(viz_tools) \
	"output/figures/FigS4.pdf" 
```

![./output/figures/FigS4.pdf](./output/figures/FigS4.pdf)

### Supplementary Figure 5
```bash
Rscript src/viz/contour_sims_plots.R \
	"output/sims/delayed/delayed_contours_simulations.csv" \
	"output/sims/immediate/immediate_contours_simulations.csv" \
	"output/sims/initial_conditions/initial_conditions_base.csv" \
	$(viz_tools) \
	"diffiso" \
	"FALSE" \
	"output/figures/FigS5.pdf"
```

![./output/figures/FigS5.pdf](./output/figures/FigS5.pdf)


### Supplementary Figure 6
```bash
Rscript src/viz/contour_sims_plots.R \
	"output/sims/delayed/delayed_contours_simulations.csv" \
	"output/sims/immediate/immediate_contours_simulations.csv" \
	"output/sims/initial_conditions/initial_conditions_base.csv" \
	$(viz_tools) \
	"both" \
	"FALSE" \
	"output/figures/FigS6.pdf"
```

![./output/figures/FigS6.pdf](./output/figures/FigS6.pdf)
