""" Log Normal noise model for transmission, applied to different
subsets of the Pakistan data. The model is fit at the national level. """

## Put the path to the riskmap lib
## in the system path
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
from riskmap3.data_process.vis_tools import *

## For measuring model performance
from sklearn.metrics import r2_score

## Some font overwrites
plt.rcParams["font.size"] = 22.

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

def low_mid_high(samples):
	low = np.percentile(samples,2.5,axis=0)
	mid = np.mean(samples,axis=0)#np.percentile(samples,50.,axis=0)
	high = np.percentile(samples,97.5,axis=0)
	return low,mid,high

if __name__ == "__main__":

	## Get the shape file
	shp = "..\\_data\\Shapefile\\Pakistan.shp"
	sf = shapefile.Reader(shp)

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
	## all states. This also drops times where we don't have case data).
	combined = pd.concat([ri,mcv2,pop,cases,br,sia],axis=1).dropna()

	## Choose relevant times (cutting off the early stuff due to time dependent
	## changes in reporting).
	combined = combined.loc(axis=0)[:,pd.to_datetime("01-01-2012"):]

	## Compute the adjusted births just to see how
	## they compare to reporting rates and cases.
	combined["adj_births"] = adj_births(combined)

	## Description
	print("\nSummary of the data and expectations for reporting rates:")
	print(combined.groupby(level=0).apply(describe))

	## Aggregate up
	## National level model
	national = combined.groupby(level=1).apply(up_sample)

	## Create model
	nat_model = TSIR(national)
	nat_model.mle(detrended=True,weighted=True,verbose=False)
	print("Inferred reporting rate = %.5f" % nat_model.reporting_rate)

	## Seasonality
	nat_model.transmission_regression(periodicity=24)
	print(nat_model.t_params)
	monthly = (nat_model.t_beta[::2] + nat_model.t_beta[1::2])/2.
	by_month = national["cases"].groupby(lambda t: t.month).sum()

	## Seasonality plot
	fig, axes = plt.subplots(figsize=(12,6))
	axes.bar(by_month.index,by_month,width=0.75,color="grey",alpha=0.75,label="Cases by Month")
	axes.set_ylabel("Cases by month",color="grey")
	axes.set_xlabel("Time of year (month)")
	axes2 = axes.twinx()
	axes2.plot(np.arange(1,13), monthly, c="C3",marker="s",ls="None",label="Inferred Seasonality")
	axes2.set_ylabel("Inferred seasonality",color="C3")
	fig.tight_layout()

	## Alpha
	alpha = nat_model.alpha
	alpha_std = np.sqrt(np.diag(nat_model.t_var)[nat_model.periodicity])
	print("Alpha estimate is {} +/- {}".format(alpha,alpha_std))

	########################## Full Forward/One step projection
	## Allocate space based on the number of sample trajectories
	## we want.
	std_logE = nat_model.std_logE
	num_samples = 10000
	full_samples = np.zeros((num_samples,nat_model.n_steps))
	full_samples_S = np.zeros((num_samples,nat_model.n_steps))
	one_step_samples = np.zeros((num_samples,nat_model.n_steps))

	## Set up initial conditions
	I_inferred = nat_model.beta[0]*(nat_model.cases+1.)-1.
	full_samples_S[:,0] = nat_model.S_bar + nat_model.Z[0]
	full_samples[:,0] = I_inferred[0]
	one_step_samples[:,0] = I_inferred[0]

	## Loop through time
	for i in range(1,len(nat_model.cases)):
		
		## Time of year for seasonality
		time_in_period = (i % nat_model.periodicity)

		## Update one step ahead and full projection lambdas
		lam = nat_model.t_beta[time_in_period]*(full_samples_S[:,i-1])*(full_samples[:,i-1]**nat_model.alpha)
		lam_one_step = nat_model.t_beta[time_in_period]*(nat_model.S_bar + nat_model.Z[i-1])*(I_inferred[i-1]**nat_model.alpha)

		## Sample for new infecteds in both cases
		I_ts = lam*np.exp(std_logE*np.random.normal(size=(num_samples,)))
		S_ts = (full_samples_S[:,i-1]+nat_model.adj_births[i-1]-I_ts)*(1. - nat_model.sia[i-1])

		## Update predictions and residuals
		full_samples_S[:,i] = S_ts
		full_samples[:,i] = I_ts
		one_step_samples[:,i] = lam_one_step*(1. + std_logE*np.random.normal(size=(num_samples,)))

	##### Plot the results
	method = 1

	## Plot the SIA lines
	fig, axes = plt.subplots(2,1,sharex=True,figsize=(13,9))
	for x in nat_model.df[nat_model.df["sia"] != 0.].index:
		axes[0].axvline(x,c="k",alpha=0.4,ls="dashed")

	## Plot either individual trajectoris (method = 0) or
	## a cloud (method = 1)
	if method == 0:
		for i in np.random.choice(num_samples,replace=False,size=(500,)):
			S = full_samples_S[i,:]
			predicted_infecteds = full_samples[i,:]
			one_step = one_step_samples[i,:]
			axes[0].plot(nat_model.df.index,S,c="C4",alpha=0.01)
			axes[1].plot(nat_model.df.index,predicted_infecteds,c="C1",alpha=0.01)
			axes[1].plot(nat_model.df.index,one_step,c="C0",alpha=0.01)
		axes[0].plot([],label="Projected S",c="C4")
		axes[1].plot([],label="Projected I",c="C1")
		axes[1].plot([],label="One-step I",c="C0")
	elif method == 1:
		full_low, full_mid, full_high = low_mid_high(full_samples)
		full_low_S, full_mid_S, full_high_S = low_mid_high(full_samples_S)
		one_step_low, one_step_mid, one_step_high = low_mid_high(one_step_samples)
		axes[0].fill_between(nat_model.df.index,full_low_S,full_high_S,color="C0",alpha=0.2)
		axes[0].plot(nat_model.df.index,full_mid_S,color="C0",label=r"S$_{t}\,|\,$S$_{0}$")
		axes[1].fill_between(nat_model.df.index,full_low,full_high,color="k",alpha=0.2)
		axes[1].fill_between(nat_model.df.index,one_step_low,one_step_high,color="C3",alpha=0.3)
		axes[1].plot(nat_model.df.index,full_mid,color="k",label=r"I$_{t}\,|\,$I$_{0}$")
		axes[1].plot(nat_model.df.index,one_step_mid,color="C3",label=r"I$_{t}\,|\,$I$_{t-1}$")
	axes[1].plot(nat_model.df.index,I_inferred,color="k",marker=".",ls="None",label="Data")
	axes[0].plot([],c="k",alpha=0.5,ls="dashed",label="SIA")
	axes[0].legend(loc=2)
	axes[1].legend(loc=2)
	axes[1].set(ylabel="Infecteds")
	axes[0].set(ylabel="Susceptibles")
	axes[0].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	axes[1].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	plt.tight_layout()
	plt.savefig("..\\_plots\\calibration_updated.png")
	

	## Print some metrics
	print(r2_score(I_inferred,one_step_mid))
	print(r2_score(I_inferred,full_mid))





	plt.show()