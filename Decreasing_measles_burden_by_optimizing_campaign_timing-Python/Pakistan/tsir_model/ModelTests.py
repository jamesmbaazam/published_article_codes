""" ModelTests.py

Test the model in three ways, one step ahead prediction, long term prediction, and 
out of sample prediction. This is simply a combined script with the tests from ex_LogNormal.py
and LogNormal.py. """
import sys
sys.path.insert(0,"..\\..\\")

## Standard imports 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shapefile

## Risk map tools
from riskmap3.tsir import TSIR
from riskmap3.map_maker import *

## For quantifying performance
from sklearn.metrics import r2_score

## Some font overwrites
plt.rcParams["font.size"] = 24.

## Some simple helper functions
#################################################################################
def adj_births(df):
	births = df.birth_rate*df.population/1000.
	try:
		mcv2 = df.mcv2
	except:
		mcv2 = 0.
	return births*(1.-0.9*df.mcv1*(1.-mcv2)-0.99*df.mcv1*mcv2)

def describe(x):
	adj_births = x.adj_births.mean()
	cases = x.cases.mean()
	rr_estimate = cases/adj_births
	return pd.Series([adj_births, cases, rr_estimate],index=["MCV Adj Births","Cases","RR Estimate"])

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

def low_mid_high(samples):
	low = np.percentile(samples,2.5,axis=0)
	mid = np.mean(samples,axis=0)#np.percentile(samples,50.,axis=0)
	high = np.percentile(samples,97.5,axis=0)
	return low,mid,high

if __name__ == "__main__":

	## Get the data
	pickle_jar = "..\\pickle_jar\\"
	ri = pd.read_pickle(pickle_jar+"extrapolated_ri.pkl").rename("mcv1")
	pop = pd.read_pickle(pickle_jar+"extrapolated_population.pkl").rename("population")
	cases = pd.read_pickle(pickle_jar+"smoothed_cases.pkl").rename("cases")
	br_wp = pd.read_pickle(pickle_jar+"extrapolated_birth_rate.pkl").rename("world_pop")
	br_dhs = pd.read_pickle(pickle_jar+"extrapolated_dhs_birth_rate.pkl").rename("dhs")
	sia = pd.read_pickle(pickle_jar+"extrapolated_sia_nofuture.pkl").rename("sia")
	mcv2 = pd.read_pickle(pickle_jar+"extrapolated_mcv2.pkl").rename("mcv2")

	## For later comparison
	updated_2017 = pd.read_pickle(pickle_jar+"updated_smoothed_cases.pkl")
	updated_cases = updated_2017.groupby("time").sum()

	## Update the cases?
	cases = updated_2017.reindex(cases.index).fillna(cases)

	## Combine birth rates estimates to take the world pop
	## where it's good and dhs otherwise.
	br_wp.loc["asia:pakistan:gilgit baltistan"] = np.nan
	br = br_wp.fillna(br_dhs).rename("birth_rate")

	## Create full DF (where we drop NA because we don't have DHS clusters in
	## all states. Could instead use interpolation here).
	## Now we also drop NaNs because this is gonna be used for model calibration,
	## and we want to pass the TSIR class only data with corresponding cases.
	combined = pd.concat([ri,mcv2,pop,cases,br,sia],axis=1).dropna()
	extrapolation = pd.concat([ri,mcv2,pop,br,sia],axis=1).dropna()

	## Choose relevant times (cutting off the early stuff due to time dependent
	## changes in reporting).
	combined = combined.loc(axis=0)[:,pd.to_datetime("01-01-2012"):pd.to_datetime("03-31-2017")]
	extrapolation = extrapolation.loc(axis=0)[:,pd.to_datetime("01-01-2012"):pd.to_datetime("01-01-2018")]

	## Compute adjusted births and nationalized dfs
	combined["adj_births"] = adj_births(combined)
	national = combined.groupby(level=1).apply(up_sample)
	n_extrap = extrapolation.groupby(level=1).apply(ex_up_sample)
	n_extrap["adj_births"] = adj_births(n_extrap)

	## Create model
	nat_model = TSIR(national)
	nat_model.mle(detrended=True,weighted=True,verbose=False)
	nat_model.transmission_regression(periodicity=24)

	########################## Extrapolation tests
	std_logE = nat_model.std_logE
	num_samples = 30000
	longterm_samples = np.zeros((num_samples,len(n_extrap)))
	longterm_samples_S = np.zeros((num_samples,len(n_extrap)))
	onestep_samples = np.zeros((num_samples,len(n_extrap)))
	onestep_samples_S = np.zeros((num_samples,len(n_extrap)))

	## Set up ICs
	I_inferred = nat_model.beta[0]*(nat_model.cases + 1.) - 1.
	longterm_samples_S[:,0] = nat_model.S_bar + nat_model.Z[0]
	longterm_samples[:,0] = I_inferred[0]
	onestep_samples_S[:,0] = nat_model.S_bar + nat_model.Z[0]
	onestep_samples[:,0] = I_inferred[0]

	## And the full infections trace
	full_cases = cases.groupby("time").sum().loc[pd.to_datetime("01-01-2012"):]
	full_I_inferred = nat_model.beta[0]*(full_cases + 1.) - 1.

	## Loop and compute samples
	for i in range(1,len(n_extrap)):
	
		## Time of year for seasonality
		time_in_period = i % nat_model.periodicity

		## Update one step ahead and full projection lambdas
		lam = nat_model.t_beta[time_in_period]*(longterm_samples_S[:,i-1])*(longterm_samples[:,i-1]**nat_model.alpha)

		## If we have data, we compute the one_step ahead projection
		if i <= nat_model.n_steps:
			lam_one_step = nat_model.t_beta[time_in_period]*(onestep_samples_S[:,i-1])*(I_inferred[i-1]**nat_model.alpha)
		else:
			lam_one_step = nat_model.t_beta[time_in_period]*(onestep_samples_S[:,i-1])*(onestep_samples[:,i-1]**nat_model.alpha)

		## Sample for the long term extrapolation
		longterm_samples[:,i] = lam*np.exp(std_logE*np.random.normal(size=(num_samples,)))
		longterm_samples_S[:,i] = (longterm_samples_S[:,i-1]+n_extrap.adj_births[i-1]-longterm_samples[:,i])*(1.-n_extrap.sia[i-1])

		## And for the step ahead (which becomes long term out of sample)
		onestep_samples[:,i] = lam_one_step*np.exp(std_logE*np.random.normal(size=(num_samples,)))
		onestep_samples_S[:,i] = (onestep_samples_S[:,i-1]+n_extrap.adj_births[i-1]-onestep_samples[:,i])*(1.-n_extrap.sia[i-1])

	## Plot the fit quality
	## Start by computing summary statistics
	lt_low, lt_mid, lt_high = low_mid_high(longterm_samples)
	lt_low_S, lt_mid_S, lt_high_S = low_mid_high(longterm_samples_S)
	os_low, os_mid, os_high = low_mid_high(onestep_samples)
	os_low_S, os_mid_S, os_high_S = low_mid_high(onestep_samples_S)

	## Initialize the figure and add SIA lines
	fig, axes = plt.subplots(2,1,sharex=True,figsize=(10,9))
	today = pd.to_datetime("13-11-2017",format="%d-%m-%Y")
	for x in n_extrap[n_extrap["sia"] != 0.].index:
		if x > today:
			color = "k"
			alpha = 0.8
		else:
			color = "grey"
			alpha = 0.3
		axes[1].axvline(x,c=color,alpha=alpha,ls="dashed")
		axes[0].axvline(x,c=color,alpha=alpha,ls="dashed")

	## Plot the model outputs
	_mode = 2
	if _mode == 1:
		axes[0].fill_between(n_extrap.index,lt_low_S,lt_high_S,color="C4",alpha=0.2)
		axes[0].plot(n_extrap.index,lt_mid_S,color="C4",label=r"S$_t$ $|$ S$_0$")
		axes[0].fill_between(n_extrap.index,os_low_S,os_high_S,color="#68829E",alpha=0.2)
		axes[0].plot(n_extrap.index,os_mid_S,color="#68829E",label=r"S$_t$ $|$ S$_{t-1}$")
		axes[1].fill_between(n_extrap.index,lt_low,lt_high,color="C1",alpha=0.2)
		axes[1].fill_between(n_extrap.index,os_low,os_high,color="#80BD9E",alpha=0.2)
		axes[1].plot(n_extrap.index,lt_mid,color="C1",label=r"I$_t$ | I$_0$")
		axes[1].plot(n_extrap.index,os_mid,color="#80BD9E",label=r"I$_t$ | I$_{t-1}$")

	elif _mode == 2:
		N = len(national)
		axes[0].fill_between(national.index,os_low_S[:N],os_high_S[:N],color="#68829E",alpha=0.2)
		axes[0].plot(national.index,os_mid_S[:N],color="#68829E",label=r"S$_t$ $|$ S$_{t-1}$")
		axes[0].fill_between(n_extrap.index[N:],os_low_S[N:],os_high_S[N:],color="#80BD9E",alpha=0.2)
		axes[0].plot(n_extrap.index[N:],os_mid_S[N:],color="#80BD9E",label=r"Out of sample S$_t$")
		axes[1].fill_between(national.index,os_low[:N],os_high[:N],color="#68829E",alpha=0.2)
		axes[1].plot(national.index,os_mid[:N],color="#68829E",label=r"I$_t$ $|$ I$_{t-1}$")
		axes[1].fill_between(n_extrap.index[N:],os_low[N:],os_high[N:],color="#80BD9E",alpha=0.2)
		axes[1].plot(n_extrap.index[N:],os_mid[N:],color="#80BD9E",label=r"Out of sample I$_t$")

	## Plot the data
	axes[1].plot(national.index,I_inferred,label="Scaled case reports",color="k",marker=".",ls="None")
	axes[1].plot(full_I_inferred.iloc[N:],label="Out of sample data",color="#FF420E",marker="s",ls="None")

	## Finish up by setting up legends
	axes[0].plot([],c="grey",alpha=0.5,ls="dashed",label="Past SIA")
	#axes[0].plot([],c="k",ls="dashed",label="Future SIA")
	axes[0].legend(loc=2)
	axes[1].legend(loc=2)
	axes[1].set(ylabel="Infected population")
	axes[0].set(ylabel="Susceptible population")
	axes[0].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	axes[1].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	axes[1].set_xlim(("03-01-2015","02-01-2018"))
	axes[1].set_xticks(["04-01-2015","08-01-2015","12-01-2015",
					 "04-01-2016","08-01-2016","12-01-2016",
					 "04-01-2017","08-01-2017","12-01-2017"])
	axes[1].set_xticklabels(["04","08\n2015","12",
						  "04","08\n2016","12",
						  "04","08\n2017","12"])


	plt.tight_layout()

	## Initialize the figure and add SIA lines
	fig, axes = plt.subplots(2,1,sharex=True,figsize=(14,9))
	today = pd.to_datetime("13-11-2017",format="%d-%m-%Y")
	for x in n_extrap[n_extrap["sia"] != 0.].index:
		if x > today:
			color = "k"
			alpha = 0.8
		else:
			color = "grey"
			alpha = 0.3
		axes[1].axvline(x,c=color,alpha=alpha,ls="dashed")
		axes[0].axvline(x,c=color,alpha=alpha,ls="dashed")

	## Plot the model outputs
	N = len(national)
	axes[0].fill_between(national.index,os_low_S[:N],os_high_S[:N],color="#68829E",alpha=0.2)
	axes[0].plot(national.index,os_mid_S[:N],color="#68829E",label=r"S$_t$ $|$ S$_{t-1}$")
	axes[0].fill_between(n_extrap.index[N:],os_low_S[N:],os_high_S[N:],color="#80BD9E",alpha=0.2)
	axes[0].plot(n_extrap.index[N:],os_mid_S[N:],color="#80BD9E",label=r"Out of sample S$_t$")
	axes[1].fill_between(national.index,os_low[:N],os_high[:N],color="#68829E",alpha=0.2)
	axes[1].plot(national.index,os_mid[:N],color="#68829E",label=r"I$_t$ $|$ I$_{t-1}$")
	axes[1].fill_between(n_extrap.index[N:],os_low[N:],os_high[N:],color="#80BD9E",alpha=0.2)
	axes[1].plot(n_extrap.index[N:],os_mid[N:],color="#80BD9E",label=r"Out of sample I$_t$")

	## Plot the data
	axes[1].plot(national.index,I_inferred,label="Scaled case reports",color="k",marker=".",ls="None")
	axes[1].plot(full_I_inferred.iloc[N:],label="Out of sample data",color="#FF420E",marker="s",ls="None")

	## Finish up by setting up legends
	axes[0].plot([],c="grey",alpha=0.5,ls="dashed",label="Past SIA")
	#axes[0].plot([],c="k",ls="dashed",label="Future SIA")
	axes[0].legend(loc=2)
	axes[1].legend(loc=2)
	axes[1].set(ylabel="Infected population")
	axes[0].set(ylabel="Susceptible population")
	axes[0].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	axes[1].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	plt.tight_layout()
	plt.savefig("..\\_plots\\out_of_sample.pdf")
	plt.show()
