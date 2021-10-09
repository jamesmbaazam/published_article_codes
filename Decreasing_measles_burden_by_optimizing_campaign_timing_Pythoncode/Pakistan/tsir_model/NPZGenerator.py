""" NPZGenerator.py

Since I keep needing to make samples of the model, this script has some tools to 
do it both at the national and at the province level. This is based on the function
from ComparePolicies2d.py """
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

def low_mid_high(samples):
	low = np.percentile(samples,2.5,axis=0)
	mid = np.mean(samples,axis=0)#np.percentile(samples,50.,axis=0)
	high = np.percentile(samples,97.5,axis=0)
	return low,mid,high

#### Extrapolation workhorses 
########################################################################################
def BaselineNPZ(pickle_jar="..\\pickle_jar\\",output_dir="_coverage_npz\\",num_samples=5000):

	""" This function does country level extrapolation with no future sia to establish a baseline. """

	## Start by getting the data that doesn't change.
	## Get the data
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

	## Create full DF (again based on ex_LogNormal.py)
	extrapolation = pd.concat([ri,mcv2,pop,br,sia],axis=1).dropna()
	extrapolation = extrapolation.loc(axis=0)[:,pd.to_datetime("01-01-2012"):]
	n_extrap = extrapolation.groupby(level=1).apply(ex_up_sample)
	n_extrap["adj_births"] = adj_births(n_extrap)

	## Prepare for extrapolation
	std_logE = nat_model.std_logE
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
	_npz_name = "no_future.npz"
	np.savez(output_dir+_npz_name,I_samples=one_step_samples,S_samples=one_step_S)

	return 

def NationalNPZs(sia_pickles,sia_coverages=np.linspace(0.2,0.6,20),
	pickle_jar="..\\pickle_jar\\",output_dir="_coverage_npz\\",num_samples=5000):

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

def ProvinceNPZs(province,sia_pickles,sia_coverages=np.linspace(0.2,0.6,20),
	pickle_jar="..\\pickle_jar\\",output_dir="_coverage_npz\\",num_samples=5000):

	""" This function does the extrapolation computations for a bunch of SIA conditions in
	serial and saves the traces to output_dir. The code is essentially copied from the function above,
	but adjusted to the province extrapolations."""

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

	## Initialize and fit the national TSIR model
	nat_model = TSIR(national)
	nat_model.mle(detrended=True,weighted=True,verbose=False)
	nat_model.transmission_regression(periodicity=24)

	## Get province level reporting rates, etc.
	## Reconstruct susceptibles, etc.
	province_name = province[province.rfind(":")+1:].replace(" ","")
	p_model = TSIR(combined.loc[province])
	p_model.mle(detrended=True,weighted=True,verbose=False)
	pop_fraction = p_model.df.population.mean()/nat_model.df.population.mean()

	## With the model fit, loop over extrapolations for
	## specific policies.
	print("\nPerforming analysis for "+province)
	for sia_pickle in sia_pickles:

		print("\nStarting computation for "+sia_pickle)
		print("...")

		for p_index, prefactor in enumerate(sia_coverage_prefactors):

			## Adjust the future SIA only to have coverage
			## set via the function inputs.
			sia = pd.read_pickle(pickle_jar+sia_pickle).rename("sia")
			sia.loc[:,"2017-12-31":] = prefactor*sia.loc[:,"2017-12-31":]

			## Create full DF (based on NatModelAtProvince.py)
			extrapolation = pd.concat([ri,mcv2,pop,br,sia],axis=1).dropna()
			extrapolation = extrapolation.loc(axis=0)[:,pd.to_datetime("01-01-2012"):]
			extrap = extrapolation.loc[province]
			extrap["adj_births"] = adj_births(extrap)

			## Prepare for extrapolation
			std_logE = nat_model.std_logE
			one_step_samples = np.zeros((num_samples,len(extrap)))
			one_step_S = np.zeros((num_samples,len(extrap)))

			## Set up ICs
			I_inferred = p_model.beta[0]*(p_model.cases + 1.) - 1.
			one_step_samples[:,0] = I_inferred[0]
			one_step_S[:,0] = pop_fraction*nat_model.S_bar + p_model.Z[0]

			## Loop through time
			for i in range(1,len(extrap)):
				
				## Time of year for seasonality
				time_in_period = i % nat_model.periodicity

				## If we have data, we compute the one_step ahead projection
				if i <= p_model.n_steps:
					lam_one_step = (nat_model.t_beta[time_in_period]/pop_fraction)*(one_step_S[:,i-1])*(I_inferred[i-1]**nat_model.alpha)
				else:
					lam_one_step = (nat_model.t_beta[time_in_period]/pop_fraction)*(one_step_S[:,i-1])*(one_step_samples[:,i-1]**nat_model.alpha)

				## Update predictions and residuals
				one_step_samples[:,i] = lam_one_step*np.exp(std_logE*np.random.normal(size=(num_samples,)))
				one_step_S[:,i] = (one_step_S[:,i-1]+extrap.adj_births[i-1]-one_step_samples[:,i])*(1. - extrap.sia[i-1])

			## Finally, save the NPZ
			_npz_name = province_name+"_"+sia_pickle[-9:-4]+"_"+str(p_index)+".npz"
			np.savez(output_dir+_npz_name,I_samples=one_step_samples,S_samples=one_step_S)

	## Create a look up table mapping p_index to 
	## a coverage value
	index_to_coverage = {i:p for i,p in enumerate(sia_coverages)}
	pickle.dump(index_to_coverage,open(output_dir+"index_to_coverage_dict.pkl","wb"))

	return 

#### Visualization (for debugging mostly)
########################################################################################
def PlotNPZ(npz_name,pickle_jar="..\\pickle_jar\\",output_dir="_coverage_npz\\"):

	""" Function to get the NPZ required and plot the susceptible and infecteds 
	distributions. Done simply, since this isn't intended to be publication worthy. """

	## Get the NPZ
	npz = np.load(output_dir+npz_name)
	infecteds = npz["I_samples"]
	susceptibles = npz["S_samples"]

	## Parse the npz name to get the full
	## SIA calendar.
	elements = npz_name.split("_")
	if npz_name[0].isdigit():
		sia_date = elements[0]+"_"+elements[1]
	else:
		sia_date = elements[1]+"_"+elements[2]
	#sia = pd.read_pickle(pickle_jar+"extrapolated_sia_"+sia_date+".pkl").rename("sia")

	## Construct the province name (handling the two exceptions).
	elements[0] = elements[0].replace("khyberpakhtoon","khyber pakhtoon")
	elements[0] = elements[0].replace("gilgitbaltistan","gilgit baltistan")

	## Subset the SIA calendar to the relevant parts
	if npz_name[0].isdigit():
		sia = sia.groupby(level=1).mean()
	else:
		sia = sia.loc["asia:pakistan:"+elements[0]]
	sia = sia.loc[pd.to_datetime("01-01-2012"):]

	## And get the case data, doing the same.
	cases = pd.read_pickle(pickle_jar+"smoothed_cases.pkl").rename("cases")
	updated_2017 = pd.read_pickle(pickle_jar+"updated_smoothed_cases.pkl")
	cases = updated_2017.reindex(cases.index).fillna(cases)
	if npz_name[0].isdigit():
		cases = cases.groupby(level=1).sum()
	else:
		cases = cases.loc["asia:pakistan:"+elements[0]]
	cases = cases.loc[pd.to_datetime("01-01-2012"):]

	## Initialize the figure and axes
	fig, axes = plt.subplots(2,1,sharex=True,figsize=(14,10))

	## Add the SIA lines 
	for x in sia[sia != 0.].index:
		axes[0].axvline(x,c="k",alpha=0.4,ls="dashed")

	## Compute the percentiles and mean
	I_low, I_mid, I_high = low_mid_high(infecteds)
	S_low, S_mid, S_high = low_mid_high(susceptibles)

	## Plot the traces
	axes[0].fill_between(sia.index,S_low,S_high,color="C0",alpha=0.2)
	axes[0].plot(sia.index,S_mid,color="C0",label=r"S$_t$ | C$_{t-1}$")
	axes[1].fill_between(sia.index,I_low,I_high,color="C4",alpha=0.2)
	axes[1].plot(sia.index,I_mid,color="C4",label=r"I$_t$ | C$_{t-1}$")

	## Plot the case data (with ball-park reporting rate, avoiding
	## re-running a full regression since this is just for debugging 
	## purposes anyway).
	rep_rate = ((I_mid[:len(cases)] + 1.)/(cases.values + 1)).mean()
	axes[1].plot(cases.index,rep_rate*(cases+1.) - 1.,label=r"Scaled C$_t$",color="k",marker=".",ls="None")

	## Legend
	axes[0].plot([],c="k",alpha=0.5,ls="dashed",label="SIA")
	axes[0].legend(loc=1)
	axes[1].legend(loc=1)
	axes[1].set(ylabel="Infecteds")
	axes[0].set(ylabel="Susceptibles")
	axes[0].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	axes[1].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	plt.tight_layout()

	return fig, axes

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
	sia_coverages = np.array([0.4]) 
	#sia_coverages = np.linspace(0.2,0.6,3)

	## And an output_dir (for where to put npz's or
	## where to get them).
	output_dir = "_updated_npz2\\"

	## Compute...
	#NationalNPZs(sia_pickles,sia_coverages,output_dir=output_dir,num_samples=20000)
	BaselineNPZ(output_dir=output_dir,num_samples=20000)
	#for p in ["punjab","sindh","islamabad","gilgit baltistan","balochistan","khyber pakhtoon"]:
	#	dn = "asia:pakistan:"+p
	#	ProvinceNPZs(dn,sia_pickles,sia_coverages,output_dir=output_dir,num_samples=10000)

	#npz_name = "no_future.npz"
	#fig, axes = PlotNPZ(npz_name,output_dir=output_dir)
	#plt.show()
