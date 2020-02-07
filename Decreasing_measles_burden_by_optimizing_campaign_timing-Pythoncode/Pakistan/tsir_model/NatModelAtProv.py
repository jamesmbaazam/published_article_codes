""" NatModelAtProv.py

Using the national level model to do province level
extrapolation."""
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

## Overwrite some risk map defaults
plt.rcParams["font.size"] = 22.

## Get the data
pickle_jar = "..\\pickle_jar\\"
ri = pd.read_pickle(pickle_jar+"extrapolated_ri.pkl").rename("mcv1")
pop = pd.read_pickle(pickle_jar+"extrapolated_population.pkl").rename("population")
cases = pd.read_pickle(pickle_jar+"smoothed_cases.pkl").rename("cases")
rejected = pd.read_pickle(pickle_jar+"rejected.pkl").rename("rejected")
br_wp = pd.read_pickle(pickle_jar+"extrapolated_birth_rate.pkl").rename("world_pop")
br_dhs = pd.read_pickle(pickle_jar+"extrapolated_dhs_birth_rate.pkl").rename("dhs")
sia = pd.read_pickle(pickle_jar+"extrapolated_sia_nofuture.pkl").rename("sia")
mcv2 = pd.read_pickle(pickle_jar+"extrapolated_mcv2.pkl").rename("mcv2")

## Save the computed traces?
province = "asia:pakistan:islamabad"
_npz_name = None #"islamabad_12_19"

## For later comparison
updated_2017 = pd.read_pickle(pickle_jar+"updated_smoothed_cases.pkl")
updated_cases = updated_2017.loc[province]

## Update the cases we use for fitting?
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
combined = combined.loc(axis=0)[:,pd.to_datetime("01-01-2012"):]#pd.to_datetime("09-24-2017")]
extrapolation = extrapolation.loc(axis=0)[:,pd.to_datetime("01-01-2012"):]

## Adjusted_births
def adj_births(df):
	births = df.birth_rate*df.population/1000.
	try:
		mcv2 = df.mcv2
	except:
		mcv2 = 0.
	return births*(1.-0.9*df.mcv1*(1.-mcv2)-0.99*df.mcv1*mcv2)
combined["adj_births"] = adj_births(combined)

## Aggregate up
## National level model
def up_sample(x):
	total_pop = x.population.sum()
	mcv1 = np.sum(x.population*x.mcv1)/total_pop
	mcv2 = np.sum(x.population*x.mcv2)/total_pop
	br = np.sum(x.population*x.birth_rate)/total_pop
	sia = np.sum(x.population*x.sia)/total_pop
	series = pd.Series([mcv1,mcv2,total_pop,x.cases.sum(),br,sia],
					  index=["mcv1","mcv2","population","cases","birth_rate","sia"])
	return series

## Create grouped and subset data frames
prov_name = province[province.rfind(":")+1:]
national = combined.groupby(level=1).apply(up_sample)
extrap = extrapolation.loc[province]
extrap["adj_births"] = adj_births(extrap)

## Create national model
nat_model = TSIR(national)
nat_model.mle(detrended=True,weighted=True,verbose=False)
nat_model.transmission_regression(periodicity=24)

## Get province level reporting rates, etc.
## Reconstruct susceptibles, etc.
p_model = TSIR(combined.loc[province])
p_model.mle(detrended=True,weighted=True,verbose=False)

########################## Full Forward/One step projection
std_logE = nat_model.std_logE
num_samples = 10000
full_samples = np.zeros((num_samples,len(extrap)))
full_samples_S = np.zeros((num_samples,len(extrap)))
one_step_samples = np.zeros((num_samples,len(extrap)))
one_step_S = np.zeros((num_samples,len(extrap)))
pop_fraction = p_model.df.population.mean()/nat_model.df.population.mean()

## Set up ICs
I_inferred = p_model.beta[0]*(p_model.cases+1.)-1.
updated_I_inferred = p_model.beta[0]*(updated_cases + 1.) - 1.
full_samples_S[:,0] = pop_fraction*nat_model.S_bar + p_model.Z[0]
full_samples[:,0] = I_inferred[0]
one_step_samples[:,0] = I_inferred[0]
one_step_S[:,0] = pop_fraction*nat_model.S_bar + p_model.Z[0]

## Loop through time
for i in range(1,len(extrap)):
	
	## Time of year for seasonality
	time_in_period = i % nat_model.periodicity

	## Update one step ahead and full projection lambdas
	lam = (nat_model.t_beta[time_in_period]/pop_fraction)*(full_samples_S[:,i-1])*(full_samples[:,i-1]**nat_model.alpha)

	## If we have data, we compute the one_step ahead projection
	if i <= p_model.n_steps:
		lam_one_step = (nat_model.t_beta[time_in_period]/pop_fraction)*(one_step_S[:,i-1])*(I_inferred[i-1]**nat_model.alpha)
	else:
		lam_one_step = (nat_model.t_beta[time_in_period]/pop_fraction)*(one_step_S[:,i-1])*(one_step_samples[:,i-1]**nat_model.alpha)

	## Sample for new infecteds in both cases
	I_ts = lam*np.exp(std_logE*np.random.normal(size=(num_samples,)))
	S_ts = (full_samples_S[:,i-1]+extrap.adj_births[i-1]-I_ts)*(1. - extrap.sia[i-1])

	## Update predictions and residuals
	full_samples_S[:,i] = S_ts
	full_samples[:,i] = I_ts
	one_step_samples[:,i] = lam_one_step*np.exp(std_logE*np.random.normal(size=(num_samples,)))
	one_step_S[:,i] = (one_step_S[:,i-1]+extrap.adj_births[i-1]-one_step_samples[:,i])*(1. - extrap.sia[i-1])
	

##### plot the results
def low_mid_high(samples):
	low = np.percentile(samples,2.5,axis=0)
	mid = np.mean(samples,axis=0)#np.percentile(samples,50.,axis=0)
	high = np.percentile(samples,97.5,axis=0)
	return low,mid,high

fig, axes = plt.subplots(2,1,sharex=True,figsize=(14,10))
## Add sias to plots
today = pd.to_datetime("13-11-2017",format="%d-%m-%Y")
for x in p_model.df[p_model.df["sia"] != 0.].index:
	axes[0].axvline(x,c="k",alpha=0.4,ls="dashed")

full_low, full_mid, full_high = low_mid_high(full_samples)
full_low_S, full_mid_S, full_high_S = low_mid_high(full_samples_S)
one_step_low, one_step_mid, one_step_high = low_mid_high(one_step_samples)
os_S_low, os_S_mid, os_S_high = low_mid_high(one_step_S)

## 2012 - 2018 extrapolation
I = len(p_model.df.index)
axes[0].fill_between(p_model.df.index,full_low_S[:I],full_high_S[:I],color="C1",alpha=0.2)
axes[0].plot(p_model.df.index,full_mid_S[:I],color="C1",label=r"S$_t$ | C$_0$")
axes[1].fill_between(p_model.df.index,full_low[:I],full_high[:I],color="C3",alpha=0.2)
axes[1].plot(p_model.df.index,full_mid[:I],color="C3",label=r"I$_t$ | C$_0$")

## 2018 and beyond projectsion (one_step from 2012 to 2018)
axes[0].fill_between(extrap.index,os_S_low,os_S_high,color="C0",alpha=0.2)
axes[0].plot(extrap.index,os_S_mid,color="C0",label=r"S$_t$ | C$_{t-1}$")
axes[1].fill_between(extrap.index,one_step_low,one_step_high,color="C4",alpha=0.2)
axes[1].plot(extrap.index,one_step_mid,color="C4",label=r"I$_t$ | C$_{t-1}$")

## data
axes[1].plot(p_model.df.index,I_inferred,label=r"Scaled C$_t$",color="k",marker=".",ls="None")
#axes[1].plot(updated_I_inferred,label="Updated data",color="C3",marker="x",ls="None")

## Legend
axes[0].plot([],c="k",alpha=0.5,ls="dashed",label="SIA")
axes[0].legend(loc=1)
axes[1].legend(loc=1)
axes[1].set(ylabel="Infecteds")
axes[0].set(ylabel="Susceptibles")
axes[0].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
axes[1].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
plt.tight_layout()
plt.savefig("..\\_plots\\"+prov_name.replace(" ","_")+"_nosia.pdf")

## Save samples
if _npz_name is not None:
	np.savez("..\\npz_jar\\"+_npz_name+"_sia.npz",I_samples=one_step_samples,S_samples=one_step_S)


## Sum the cases after a given time
cut_off_time = pd.to_datetime("31-12-2017",format="%d-%m-%Y")
I = np.argmin(np.abs(extrap.index - cut_off_time))
future_cases = np.sum(one_step_mid[I:])
print("Infections after 2018 in this model = %.0f" % future_cases)
plt.show()