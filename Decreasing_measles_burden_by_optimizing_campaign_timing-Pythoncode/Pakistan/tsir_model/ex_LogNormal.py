""" Extrapolated version of the log normal script. Extrapolated as in
we let the model continue past where there's data."""
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
plt.rcParams["font.size"] = 26.

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

## For later comparison
updated_2017 = pd.read_pickle(pickle_jar+"updated_smoothed_cases.pkl")
updated_cases = updated_2017.groupby("time").sum()

## Update the cases?
cases = updated_2017.reindex(cases.index).fillna(cases)

## Save the computed traces?
_npz_name = None #"12_19"

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
def ex_up_sample(x):
	total_pop = x.population.sum()
	mcv1 = np.sum(x.population*x.mcv1)/total_pop
	mcv2 = np.sum(x.population*x.mcv2)/total_pop
	br = np.sum(x.population*x.birth_rate)/total_pop
	sia = np.sum(x.population*x.sia)/total_pop
	series = pd.Series([mcv1,mcv2,total_pop,br,sia],
					  index=["mcv1","mcv2","population","birth_rate","sia"])
	return series
national = combined.groupby(level=1).apply(up_sample)
n_extrap = extrapolation.groupby(level=1).apply(ex_up_sample)
n_extrap["adj_births"] = adj_births(n_extrap)

## Create model
nat_model = TSIR(national)
nat_model.mle(detrended=True,weighted=True,verbose=False)
sig2_rho = np.diag(nat_model.beta_var)[0]
sig2_p = sig2_rho/(nat_model.beta[0]**4)#  - 0.25*(sig2_rho**2)/(nat_model.beta[0]**6)
print("Inferred reporting rate = {} +/- {}".format(nat_model.reporting_rate,np.sqrt(sig2_p)))

## Seasonality
#### Seasonality
nat_model.transmission_regression(periodicity=24)
monthly = (nat_model.t_beta[::2] + nat_model.t_beta[1::2])/2.
by_month = national["cases"].groupby(lambda t: t.month).sum()
by_month_2013 = national.loc["01-01-2013":"31-12-2013","cases"].groupby(lambda t: t.month).sum()
by_month_2017 = national.loc["01-01-2017":"31-12-2017","cases"].groupby(lambda t: t.month).sum()
periodic_monthly = [monthly[-1]] + monthly.tolist() + [monthly[0]]

## Seasonality uncertainty
## Gotten via
## beta = f(x,Sbar) + df/dx*sig_x + df/dSbar*sig_Sbar
sig2s = np.diag(nat_model.t_var)
sig2 = sig2s[:nat_model.periodicity] + sig2s[nat_model.periodicity+1]/(nat_model.S_bar**2)
sig = np.exp(nat_model.t_params[:nat_model.periodicity])*np.sqrt(sig2)/nat_model.S_bar
low = nat_model.t_beta
high = nat_model.t_beta
low = (low[::2] + low[1::2])/2.
high = (high[::2] + high[1::2])/2.
sig = np.sqrt((sig[::2]**2 + sig[1::2]**2)/2.)
low = low - sig
high = high + sig
periodic_low = [low[-1]] + low.tolist() + [low[0]]
periodic_high = [high[-1]] + high.tolist() + [high[0]]

## Seasonality plot
method = 1
fig, axes = plt.subplots(figsize=(9,11))

## Split bars
if method == 0:
	axes.bar(np.arange(1,13)+0.75/4,by_month,width=0.75/2,color="grey",alpha=0.75)#,label="2012-18 cases")
	axes.bar(np.arange(1,13)-0.75/4,by_month_2013.as_matrix()+by_month_2017.as_matrix(),width=0.75/2,color="C1",alpha=0.7,label="2017 cases")
	axes.bar(np.arange(1,13)-0.75/4,by_month_2013,width=0.75/2,color="C0",alpha=0.9,label="2013 cases")
	axes.set_ylabel("Cases by month",color="grey")
	axes.set_ylim((0.,3000.))
	axes.set_xlabel("Time of year (month)")
	axes.set_xlim((0.5,12.5))
	axes2 = axes.twinx()
	#axes2.plot(np.arange(1,13), monthly, c="C3",ls="steps",lw=3,label="Inferred Seasonality")
	axes2.plot(np.arange(0,14), periodic_monthly, c="k",ls="steps",lw=3,label="Inferred Seasonality")
	axes2.set_ylabel("Inferred seasonality",color="k")
	axes2.set_xlim((0.5,12.5))
	axes.legend()
	fig.tight_layout()

## Layered bars
elif method == 1:
	#axes.grid(color="grey",alpha=0.2)
	axes.bar(np.arange(1,13),by_month,width=0.75,color="grey",alpha=1,label="Cases by month")
	#axes.bar(np.arange(1,13),by_month_2013.as_matrix()+by_month_2017.as_matrix(),width=0.75,color="C1",alpha=0.7,label="2017 cases")
	#axes.bar(np.arange(1,13),by_month_2013,width=0.75,color="C0",alpha=0.9,label="2013 cases")
	axes.set_ylabel("Cases by month",color="k")
	axes.set_ylim((0.,3000.))
	axes.set_xlabel("Time of year (month)")
	axes.set_xlim((0.5,12.5))
	#axes.text(10,2750,"Cases by month",color="grey",horizontalalignment="left",verticalalignment="center",fontsize=20,
	#		bbox={"facecolor":"white","alpha":0.75,"pad":7,"edgecolor":"None"})
	axes2 = axes.twinx()
	axes2.fill_between(np.arange(0,14),periodic_low,periodic_high,color="C3",alpha=0.2,step="post")
	axes2.plot(np.arange(0,14), periodic_monthly, c="C3",lw=3,label="Inferred seasonality",ls="steps-post")
	#axes2.plot(np.arange(0,14),[2e-7]*14,lw=5,c="C1")
	axes2.plot(np.arange(5,11),[1.75e-7]*6,lw=6,c="C0")
	axes2.set_ylabel(r"Monthly $\beta_t$",color="k")
	axes2.set_xlim((0.5,12.5))
	axes2.set_ylim((1.5e-7,7.5e-7))
	axes2.ticklabel_format(axis="y",style="sci",scilimits=(0,1))
	axes.plot([],c="C3",lw=3,label="Inferred seasonality")
	axes.plot([],c="C0",lw=6,label="Inferred low season")
	axes.set_xticks(np.arange(2,13,2))
	#axes.set_xticks(np.arange(1,13))
	#axes.set_xticklabels(["J","F","M","A","M","J","J","A","S","O","N","D"])
	axes.legend()
	fig.tight_layout()
	print(monthly)

plt.savefig("..\\_plots\\national_seasonality_v6.png")


########################## Full Forward/One step projection
std_logE = nat_model.std_logE
num_samples = 10000
full_samples = np.zeros((num_samples,len(n_extrap)))
full_samples_S = np.zeros((num_samples,len(n_extrap)))
one_step_samples = np.zeros((num_samples,len(n_extrap)))
one_step_S = np.zeros((num_samples,len(n_extrap)))

## Set up ICs
I_inferred = nat_model.beta[0]*(nat_model.cases + 1.) - 1.
updated_I_inferred = nat_model.beta[0]*(updated_cases + 1.) - 1.
full_samples_S[:,0] = nat_model.S_bar + nat_model.Z[0]
full_samples[:,0] = I_inferred[0]
one_step_samples[:,0] = I_inferred[0]
one_step_S[:,0] = nat_model.S_bar + nat_model.Z[0]

## Loop through time
for i in range(1,len(n_extrap)):
	
	## Time of year for seasonality
	time_in_period = i % nat_model.periodicity

	## Update one step ahead and full projection lambdas
	lam = nat_model.t_beta[time_in_period]*(full_samples_S[:,i-1])*(full_samples[:,i-1]**nat_model.alpha)

	## If we have data, we compute the one_step ahead projection
	if i <= nat_model.n_steps:
		lam_one_step = nat_model.t_beta[time_in_period]*(one_step_S[:,i-1])*(I_inferred[i-1]**nat_model.alpha)
	else:
		lam_one_step = nat_model.t_beta[time_in_period]*(one_step_S[:,i-1])*(one_step_samples[:,i-1]**nat_model.alpha)

	## Sample for new infecteds in both cases
	I_ts = lam*np.exp(std_logE*np.random.normal(size=(num_samples,)))
	S_ts = (full_samples_S[:,i-1]+n_extrap.adj_births[i-1]-I_ts)*(1. - n_extrap.sia[i-1])

	## Update predictions and residuals
	full_samples_S[:,i] = S_ts
	full_samples[:,i] = I_ts
	if i <= nat_model.n_steps:
		one_step_samples[:,i] = lam_one_step*np.exp(std_logE*np.random.normal(size=(num_samples,)))
		one_step_S[:,i] = (one_step_S[:,i-1]+n_extrap.adj_births[i-1]-one_step_samples[:,i])*(1. - n_extrap.sia[i-1])
	else:
		one_step_samples[:,i] = lam_one_step*np.exp(std_logE*np.random.normal(size=(num_samples,)))
		one_step_S[:,i] = (one_step_S[:,i-1]+n_extrap.adj_births[i-1]-one_step_samples[:,i])*(1. - n_extrap.sia[i-1])

##### plot the results
method = 1
def low_mid_high(samples):
	low = np.percentile(samples,2.5,axis=0)
	mid = np.mean(samples,axis=0)#np.percentile(samples,50.,axis=0)
	high = np.percentile(samples,97.5,axis=0)
	return low,mid,high

fig, axes = plt.subplots(2,1,sharex=True,figsize=(14,9))
## Add sias to plots
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

## Plot results
if method == 0:
	for i in np.random.choice(num_samples,replace=False,size=(500,)):
		S_os = one_step_S[i]
		S = full_samples_S[i]
		predicted_infecteds = full_samples[i]
		one_step = one_step_samples[i]
		axes[0].plot(n_extrap.index,S_os,c="C3",alpha=0.01)
		axes[0].plot(n_extrap.index,S,c="C4",alpha=0.01)
		axes[1].plot(n_extrap.index,predicted_infecteds,c="C1",alpha=0.01)
		axes[1].plot(n_extrap.index,one_step,c="C0",alpha=0.01)
	axes[0].plot([],label="Projected S",c="C4")
	axes[1].plot([],label="Projected I",c="C1")
	axes[1].plot([],label="One-step I",c="C0")
elif method == 1:
	full_low, full_mid, full_high = low_mid_high(full_samples)
	full_low_S, full_mid_S, full_high_S = low_mid_high(full_samples_S)
	one_step_low, one_step_mid, one_step_high = low_mid_high(one_step_samples)
	os_S_low, os_S_mid, os_S_high = low_mid_high(one_step_S)
	#axes[0].fill_between(n_extrap.index,full_low_S,full_high_S,color="C4",alpha=0.2)
	#axes[0].plot(n_extrap.index,full_mid_S,color="C4",label="Projected S")
	axes[0].fill_between(n_extrap.index,os_S_low,os_S_high,color="#68829E",alpha=0.2)
	axes[0].plot(n_extrap.index,os_S_mid,color="#68829E",label="Projected S")
	#axes[1].fill_between(n_extrap.index,full_low,full_high,color="C1",alpha=0.2)
	axes[1].fill_between(n_extrap.index,one_step_low,one_step_high,color="#80BD9E",alpha=0.2)
	#axes[1].plot(n_extrap.index,full_mid,color="C1",label="Projected I")
	axes[1].plot(n_extrap.index,one_step_mid,color="#80BD9E",label="Projected I")
axes[1].plot(national.index,I_inferred,label="Scaled case reports",color="k",marker=".",ls="None")
axes[1].plot(updated_I_inferred,label="Out of sample data",color="#FF420E",marker="s",ls="None")
axes[0].plot([],c="grey",alpha=0.5,ls="dashed",label="Past SIA")
#axes[0].plot([],c="k",ls="dashed",label="Future SIA")
axes[0].legend(loc=2)
axes[1].legend(loc=2)
axes[1].set(ylabel="Infecteds")
axes[0].set(ylabel="Susceptibles")
axes[0].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
axes[1].ticklabel_format(axis="y",style="sci",scilimits=(0,1))
plt.tight_layout()
plt.savefig("..\\_plots\\out_of_sample_perf.png")

## Extrapolation zoom
fig, axes = plt.subplots(figsize=(12,8))
today = pd.to_datetime("13-11-2017",format="%d-%m-%Y")
for x in n_extrap[n_extrap["sia"] != 0.].index:
	if x > today:
		color = "k"
		alpha = 0.8
	else:
		color = "grey"
		alpha = 0.3
	axes.axvline(x,c=color,alpha=alpha,ls="dashed")
axes.fill_between(n_extrap.index,one_step_low,one_step_high,color="C4",alpha=0.2)
axes.plot(n_extrap.index,one_step_mid,color="C4",label="Projected I")
axes.plot(national.index,I_inferred,label="Scaled case reports",color="k",marker=".",ls="None")
axes.plot(updated_I_inferred.loc["04-01-2017":],label="Updated data",color="C3",marker="x",ls="None")
axes.legend(loc=2)
axes.set(xlabel="Time",xlim=(None,"01-01-2018"))
plt.tight_layout()

## Save samples
if _npz_name is not None:
	np.savez("..\\npz_jar\\"+_npz_name+"_sia.npz",I_samples=one_step_samples,S_samples=one_step_S)

## Sum the cases after a given time
cut_off_time = pd.to_datetime("31-12-2017",format="%d-%m-%Y")
I = np.argmin(np.abs(n_extrap.index - cut_off_time))
future_cases_samples = one_step_samples[:,I:]
total_samples = np.sum(future_cases_samples,axis=1)
future_cases = total_samples.mean()
future_cases_std = np.std(total_samples)
print("Infections after 2018 in this model = %.0f +/- %.0f" % (future_cases,future_cases_std))
plt.show()
