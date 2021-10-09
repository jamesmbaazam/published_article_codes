""" ComparePolicies2d.py

This is a set of scripts based on DynamicsSIAEst.py and ex_LogNormal.py
to create npz files for lots of policies, import them, and visualize them
in 2d coverage, timing space. """
import sys
sys.path.insert(0,"..\\..\\")

## Standard imports 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shapefile
import pickle

## Risk map tools
from riskmap3.tsir import TSIR
from riskmap3.map_maker import *
from riskmap3.data_process.vis_tools import *

## Overwrite some risk map defaults
plt.rcParams["font.size"] = 24.

#### Basic data processing functions
########################################################################################
def adj_births(df):
	births = df.birth_rate*df.population/1000.
	try:
		mcv2 = df.mcv2
	except:
		mcv2 = 0.
	return births*(1.-0.9*df.mcv1*(1.-mcv2)-0.99*df.mcv1*mcv2)

def up_sample(x):
	total_pop = x.population.sum()
	mcv1 = np.sum(x.population*x.mcv1)/total_pop
	mcv2 = np.sum(x.population*x.mcv2)/total_pop
	br = np.sum(x.population*x.birth_rate)/total_pop
	sia = np.sum(x.population*x.sia)/total_pop
	series = pd.Series([mcv1,mcv2,total_pop,x.cases.sum(),br,sia],
					  index=["mcv1","mcv2","population","cases","birth_rate","sia"])
	return series

def ex_up_sample(x):
	total_pop = x.population.sum()
	mcv1 = np.sum(x.population*x.mcv1)/total_pop
	mcv2 = np.sum(x.population*x.mcv2)/total_pop
	br = np.sum(x.population*x.birth_rate)/total_pop
	sia = np.sum(x.population*x.sia)/total_pop
	series = pd.Series([mcv1,mcv2,total_pop,br,sia],
					  index=["mcv1","mcv2","population","birth_rate","sia"])
	return series

def describe_future(samples,
					cut_off_time=pd.to_datetime("31-12-2017",format="%d-%m-%Y"),
					time_index = pd.DatetimeIndex(start="15-01-2012",end="31-12-2020",freq="SM")):
	cut_off_I = np.argmin(np.abs(time_index - cut_off_time))
	future_cases_samples = samples[:,cut_off_I:]
	total_samples = np.sum(future_cases_samples,axis=1)
	future_cases = total_samples.mean()
	future_cases_std = np.std(total_samples)
	return future_cases, future_cases_std

#### Extrapolation workhorse
########################################################################################
def MakeNPZs(sia_pickles,sia_coverages=np.linspace(0.2,0.6,20),
			 pickle_jar="..\\pickle_jar\\",output_dir="_coverage_npz\\"):

	""" This function does the extrapolation computations for a bunch of SIA conditions in
	serial and saves the traces to output_dir. The code is essentially copied from DynamicSIAEst.py"""

	## This is gonna be done is a really dumpy way, looping through
	## copy-pasted iterations of LogNormal.py with a prefactor on SIA 
	## coverage, storing the model outputs for each.
	## Base coverage is 0.4, so prefactor implies 1.1*0.4 = .44 base etc.
	sia_coverage_prefactors = sia_coverages/0.4

	## Start by getting the data that doesn't change.
	## Get the data
	pickle_jar = "..\\pickle_jar\\"
	ri = pd.read_pickle(pickle_jar+"extrapolated_ri.pkl").rename("mcv1")
	pop = pd.read_pickle(pickle_jar+"extrapolated_population.pkl").rename("population")
	cases = pd.read_pickle(pickle_jar+"smoothed_cases.pkl").rename("cases")
	br_wp = pd.read_pickle(pickle_jar+"extrapolated_birth_rate.pkl").rename("world_pop")
	br_dhs = pd.read_pickle(pickle_jar+"extrapolated_dhs_birth_rate.pkl").rename("dhs")
	mcv2 = pd.read_pickle(pickle_jar+"extrapolated_mcv2.pkl").rename("mcv2")
	br_wp.loc["asia:pakistan:gilgit baltistan"] = np.nan
	br = br_wp.fillna(br_dhs).rename("birth_rate")

	## For later comparison
	updated_2017 = pd.read_pickle(pickle_jar+"updated_smoothed_cases.pkl")
	updated_cases = updated_2017.groupby("time").sum()

	## Update the cases?
	cases = updated_2017.reindex(cases.index).fillna(cases)

	## For now, use the no future SIA to parameterize the model
	sia = pd.read_pickle(pickle_jar+"extrapolated_sia_nofuture.pkl").rename("sia")

	## Create the dataframe for fitting the model, based on the 
	## procedure in ex_LogNormal.py
	combined = pd.concat([ri,mcv2,pop,cases,br,sia],axis=1).dropna()
	combined = combined.loc(axis=0)[slice(None),pd.to_datetime("01-01-2012"):]
	combined["adj_births"] = adj_births(combined)
	national = combined.groupby(level=1).apply(up_sample)

	## Initialize and fit the TSIR model
	nat_model = TSIR(national)
	nat_model.mle(detrended=True,weighted=True,verbose=False)
	nat_model.transmission_regression(periodicity=24)

	## With the model fit, loop over extrapolations for
	## specific policies.
	for sia_pickle in sia_pickles:

		print("\nStarting computation for "+sia_pickle)
		print("...")

		for p_index, prefactor in enumerate(sia_coverage_prefactors):

			## Adjust the future SIA only to have coverage
			## set via the function inputs.
			sia = pd.read_pickle(pickle_jar+sia_pickle).rename("sia")
			sia.loc[:,"2017-12-31":] = prefactor*sia.loc[:,"2017-12-31":]

			## Create full DF (again based on ex_LogNormal.py)
			extrapolation = pd.concat([ri,mcv2,pop,br,sia],axis=1).dropna()
			extrapolation = extrapolation.loc(axis=0)[:,pd.to_datetime("01-01-2012"):]
			n_extrap = extrapolation.groupby(level=1).apply(ex_up_sample)
			n_extrap["adj_births"] = adj_births(n_extrap)

			## Prepare for extrapolation
			std_logE = nat_model.std_logE
			num_samples = 5000
			one_step_samples = np.zeros((num_samples,len(n_extrap)))
			one_step_S = np.zeros((num_samples,len(n_extrap)))

			## Set up ICs
			I_inferred = nat_model.beta[0]*(nat_model.cases + 1.) - 1.
			one_step_samples[:,0] = I_inferred[0]
			one_step_S[:,0] = nat_model.S_bar + nat_model.Z[0]

			## Loop through time
			for i in range(1,len(n_extrap)):
				
				## Time of year for seasonality
				time_in_period = i % nat_model.periodicity

				## If we have data, we compute the one_step ahead projection
				if i <= nat_model.n_steps:
					lam = nat_model.t_beta[time_in_period]*(one_step_S[:,i-1])*(I_inferred[i-1]**nat_model.alpha)
				else:
					lam = nat_model.t_beta[time_in_period]*(one_step_S[:,i-1])*(one_step_samples[:,i-1]**nat_model.alpha)

				## Update predictions and residuals
				if i <= nat_model.n_steps:
					one_step_samples[:,i] = lam*np.exp(std_logE*np.random.normal(size=(num_samples,)))
					one_step_S[:,i] = (one_step_S[:,i-1]+n_extrap.adj_births[i-1]-one_step_samples[:,i])*(1. - n_extrap.sia[i-1])
				else:
					one_step_samples[:,i] = lam*np.exp(std_logE*np.random.normal(size=(num_samples,)))
					one_step_S[:,i] = (one_step_S[:,i-1]+n_extrap.adj_births[i-1]-one_step_samples[:,i])*(1. - n_extrap.sia[i-1])

			## Finally, save the NPZ
			_npz_name = sia_pickle[-9:-4]+"_"+str(p_index)+".npz"
			np.savez(output_dir+_npz_name,I_samples=one_step_samples,S_samples=one_step_S)

	## Create a look up table mapping p_index to 
	## a coverage value
	index_to_coverage = {i:p for i,p in enumerate(sia_coverages)}
	pickle.dump(index_to_coverage,open(output_dir+"index_to_coverage_dict.pkl","wb"))

	return #one_step_samples, one_step_S, n_extrap, national, I_inferred

if __name__ == "__main__":

	## Get the list of SIA pickles
	years = ["18","19"]
	months = [str(i) for i in range(1,13)]
	sia_pickles = []
	for y in years:
		for m in months:
			pickle_name = "extrapolated_sia_"
			if len(m) == 1:
				pickle_name = pickle_name+"0"
			sia_pickles.append(pickle_name+m+"_"+y+".pkl")

	## Set up the coverages you want
	#sia_coverages = np.array([0.4]) 
	sia_coverages = np.linspace(0.2,0.6,20)

	## Pick the mode for this script
	_mode = "visualize" 

	## And an output_dir (for where to put npz's or
	## where to get them).
	output_dir = "_coverage_npz2\\"

	## If we're computing, use the NPZ function
	## to make hella npzs
	if _mode == "compute":
		MakeNPZs(sia_pickles,sia_coverages,output_dir=output_dir)
		sys.exit()

	## Load the outputs, computing the end result in each case
	## Start by allocating storage.
	coverage = []
	timing = []
	future_cases = []
	future_std = []

	## And getting the look up table
	index_to_coverage = pickle.load(open(output_dir+"index_to_coverage_dict.pkl","rb"))

	## Loop over computations
	for sia_pickle in sia_pickles:

		## Get the timing information
		date_string = sia_pickle[-9:-4]
		this_timing = int(date_string[:2])+12*(int(date_string[3:])-18)

		for i in range(len(sia_coverages)):

			## Construct the NPZ file name
			_npz_name = date_string+"_"+str(i)+".npz"

			## Load the NPZ
			npz = np.load(output_dir+_npz_name)

			## Compute the future cases
			cases, std = describe_future(npz["I_samples"])

			## Get the coverage associated with
			## this file.
			this_coverage = index_to_coverage[i]

			## Store the results
			coverage.append(this_coverage)
			timing.append(this_timing)
			future_cases.append(cases)
			future_std.append(std)

	## Convert them to np arrays
	coverage = np.array(coverage)
	timing = np.array(timing)
	future_cases = np.array(future_cases)
	future_std = np.array(future_std)

	## Reshape the results
	future_cases_m = np.reshape(future_cases,(24,len(sia_coverages)))

	## Make a heat map
	fig, axes = plt.subplots(figsize=(12,10))

	## Imshow for the matrix
	cax = axes.imshow(future_cases_m.T,
					  origin="lower",
					  extent=(1,24,sia_coverages[0],sia_coverages[-1]),
					  cmap="viridis_r",
					  aspect="auto")

	## Set up the labeling
	axes.set_xticks(np.arange(2,25,2))
	axes.set_xticklabels(["2","4","6\n2018","8","10","12"]+["2","4","6\n2019","8","10","12"])
	axes.set_title("Average total infections, 2018-21")

	## Set up y axis labels
	y_ticks = np.array([0.2,0.3,0.4,0.5,0.6])
	axes.set_yticks(y_ticks)
	axes.set_yticklabels(["{:.1f} ({:.1f})".format(y,y/0.4) for y in y_ticks])
	axes.set_ylabel("SIA coverage (value relative to 2013-14)")

	## Create a colorbar
	cbar = fig.colorbar(cax)
	cbar.formatter.set_powerlimits((0, 0))
	cbar.ax.yaxis.set_offset_position('left')   
	cbar.update_ticks()

	## Finish up
	plt.tight_layout()

	## Scatterplot version
	fig, axes = plt.subplots(figsize=(12,10))

	## Plot the low seasons
	axes.axvspan(5, 10, alpha=0.2, color="grey")
	axes.axvspan(12+5, 12+10, alpha=0.2, color="grey")

	## Plot the data
	cax = axes.scatter(timing,future_cases,marker="o",c=coverage,
				 	   cmap=plt.cm.get_cmap("viridis"),s=15**2)

	## Set up the colorbar
	cbar_ticks = np.array([0.2,0.3,0.4,0.5,0.6])
	cbar = fig.colorbar(cax,ticks=cbar_ticks)
	cbar.ax.set_yticklabels(["{:.1f} ({:.1f})".format(c,c/0.4) for c in cbar_ticks])
	cbar.ax.set_ylabel("SIA efficacy", rotation=-90, va="bottom")

	# Set up the labeling
	axes.set_xticks(np.arange(2,25,2))
	axes.set_xticklabels(["2","4","6\n2018","8","10","12"]+["2","4","6\n2019","8","10","12"])
	axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	axes.set_ylabel("Expected total infections, 2018-21")
	axes.set_xlabel("SIA timing (month)")

	## Some text for the low seasons
	axes.text(7.5, 4.575e6, "2018 low\nseason",
			verticalalignment='center', horizontalalignment='center',
			color="k", fontsize=22)
	axes.text(12.+7.5, 4.575e6, "2019 low\nseason",
			verticalalignment='center', horizontalalignment='center',
			color="k", fontsize=22)

	## Finish up
	plt.tight_layout()
	plt.savefig("..\\_plots\\2d_optim.pdf")

	## Scatterplot version (Greyed with highlights)
	fig, axes = plt.subplots(figsize=(12,10))

	## Plot the low seasons
	axes.axvspan(5, 10, alpha=0.2, color="grey")
	axes.axvspan(12+5, 12+10, alpha=0.2, color="grey")

	## Plot the data
	cax = axes.scatter(timing,future_cases,marker="o",c=coverage,
				 	   cmap=plt.cm.get_cmap("gray"),s=15**2)

	## Set up the colorbar
	cbar_ticks = np.array([0.2,0.3,0.4,0.5,0.6])
	cbar = fig.colorbar(cax,ticks=cbar_ticks)
	cbar.ax.set_yticklabels(["{:.1f} ({:.1f})".format(c,c/0.4) for c in cbar_ticks])
	cbar.ax.set_ylabel("SIA efficacy", rotation=-90, va="bottom")

	# Set up the labeling
	axes.set_xticks(np.arange(2,25,2))
	axes.set_xticklabels(["2","4","6\n2018","8","10","12"]+["2","4","6\n2019","8","10","12"])
	axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	axes.set_ylabel("Expected total infections, 2018-21")
	axes.set_xlabel("SIA timing (month)")

	## Some text for the low seasons
	axes.text(7.5, 4.575e6, "2018 low\nseason",
			verticalalignment='center', horizontalalignment='center',
			color="k", fontsize=22)
	axes.text(12.+7.5, 4.575e6, "2019 low\nseason",
			verticalalignment='center', horizontalalignment='center',
			color="k", fontsize=22)

	## Highlight particular coverages (choose indices)
	mean_cov = 9
	low_cov = 5
	high_cov = 14
	x0 = np.where(coverage == sia_coverages[low_cov])[0]
	x1 = np.where(coverage == sia_coverages[mean_cov])[0]
	x2 = np.where(coverage == sia_coverages[high_cov])[0]
	y = np.where(timing == 10)[0]
	#x = np.intersect1d(x,y)
	axes.plot(timing[x0],future_cases[x0],color="C3",lw=4,label=r"30% efficacy")
	axes.plot(timing[x1],future_cases[x1],color="C1",lw=4,label=r"40% efficacy")
	axes.plot(timing[x2],future_cases[x2],color="C2",lw=4,label=r"50% efficacy")
		
	## Make the legend
	_type = 0

	if _type == 0:
		axes.legend()

	elif _type == 1:	
		axes.text(17, 2.5e6, r"30% efficacy",
				verticalalignment='center', horizontalalignment='left',
				color="C3", fontsize=30)
		axes.text(17, 2.25e6, r"40% efficacy",
				verticalalignment='center', horizontalalignment='left',
				color="C1", fontsize=30)
		axes.text(17, 2.e6, r"50% efficacy",
				verticalalignment='center', horizontalalignment='left',
				color="C2", fontsize=30)

	## Finish up
	plt.tight_layout()
	plt.savefig("..\\_plots\\2d_optim_highlights.pdf")

	## Print a summary
	times_of_interest = [10,15,18]
	y = np.where(np.isin(timing,times_of_interest))[0]
	for cov in sia_coverages[[low_cov,mean_cov,high_cov]]:
		print("\nFor coverage = {} ...".format(cov))
		x = np.intersect1d(np.where(coverage == cov)[0],y)
		for t, f, s in zip(timing[x],future_cases[x],future_std[x]):
			print("Timing {} has {} +/- {} cases".format(t,f,s))








	plt.show()










	