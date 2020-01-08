""" LightsAndSeasonality.py

Do annualized fluctuations in major cities' night time lights correlate with
measles seasonality? """
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

## Some simple helper functions
#################################################################################
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

def GetSeasonality():

	""" Subroutine to partially fit a TSIR model and get seasonality. """

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

	## Aggregate up
	## National level model
	national = combined.groupby(level=1).apply(up_sample)

	## Create model
	nat_model = TSIR(national)
	nat_model.mle(detrended=True,weighted=True,verbose=False)

	## Seasonality
	nat_model.transmission_regression(periodicity=24)
	monthly = (nat_model.t_beta[::2] + nat_model.t_beta[1::2])/2.

	## And the uncertainty
	sig2s = np.diag(nat_model.t_var)
	sig2 = sig2s[:nat_model.periodicity] + sig2s[nat_model.periodicity+1]/(nat_model.S_bar**2)
	sig = np.exp(nat_model.t_params[:nat_model.periodicity])*np.sqrt(sig2)/nat_model.S_bar
	sig = np.sqrt((sig[::2]**2 + sig[1::2]**2)/2.)

	## Convert them to series
	monthly = pd.Series(monthly,index=np.arange(1,13),name="beta")
	sig = pd.Series(sig,index=np.arange(1,13),name="beta_sig")

	return monthly, sig

if __name__ == "__main__":

	## Get the seasonality
	beta, sig_beta = GetSeasonality()

	## Get the night time lights
	df = pd.read_pickle("..\\pickle_jar\\nighttimelights.pkl")

	## Groupby month of the year for brightness traces
	df["per_pixel"] = df["total"]/df["num_pixels"]
	df = df.groupby(["dot_name",lambda s: s[1].month]).agg({"per_pixel":["mean","var"],
															"total":"mean",
															"num_pixels":"mean"})
	summary = df["per_pixel"].groupby(level=1).mean()
	mean = summary["mean"].rename("brightness")
	std = np.sqrt(summary["var"]).rename("brightness_sig")

	## Make a comprehensive dataframe
	monthly = pd.concat([mean,std,beta,sig_beta],axis=1)
	print(monthly)
	corr = monthly[["brightness","beta"]].corr()
	print(corr)

	## Fit a linear regression
	p, V = np.polyfit(mean.values,beta.values,deg=1,cov=True)
	
	## Sample it
	p_samples = np.random.multivariate_normal(mean=p,cov=V,size=(5000,))
	x = np.linspace(0.8*mean.min(),1.15*mean.max(),100)
	samples = np.zeros((5000,100))
	for i, p_i in enumerate(p_samples):
		samples[i] = p_i[0]*x + p_i[1]
	low = np.percentile(samples,2.5,axis=0)
	high = np.percentile(samples,97.5,axis=0)
	mid = p[0]*x + p[1]

	## Make a plot
	fig, axes = plt.subplots(figsize=(12,10))
	axes.grid(color="grey",alpha=0.2)
	axes.errorbar(mean.values,beta.values,yerr=sig_beta.values,xerr=std.values,
				  ls="None",color="k",marker="o",markersize=15)
	axes.fill_between(x,low,high,color="grey",alpha=0.3)
	axes.plot(x,mid,color="k",lw=2)
	axes.set(ylabel=r"Monthly $\beta_t$",xlabel="Monthly nighttime light brightness in Karachi and Lahore")
	axes.text(16,5.8e-7,"Correlation = {0:0.3f}".format(corr.values[0,1]),fontsize=32,
			  horizontalalignment="left",verticalalignment="top",color="xkcd:red wine")
	plt.tight_layout()
	plt.savefig("..\\_plots\\transmission_vs_lights.png")


	fig, axes = plt.subplots(figsize=(12,10))
	axes.grid(color="grey",alpha=0.2)
	axes.errorbar(mean.values,beta.values,yerr=sig_beta.values,xerr=std.values,
				  ls="None",color="k",marker="o",markersize=15)
	axes.fill_between(x,low,high,color="grey",alpha=0.3)
	axes.plot(x,mid,color="k",lw=2)
	axes.set(ylabel=r"Monthly $\beta_t$",xlabel="Monthly nighttime light brightness in Karachi and Lahore")
	axes.text(35.5,2.e-7,"Correlation = {0:0.3f}".format(corr.values[0,1]),fontsize=32,
			  horizontalalignment="left",verticalalignment="bottom",color="xkcd:red wine")
	plt.tight_layout()

	## Add an inset for the timeseries
	inset_dim = [0.15, 0.685, 0.3, 0.22] ## left, bottom, width, height
	axes2 = fig.add_axes(inset_dim)

	## Make a plot of night time lights over time
	for dn, subframe in df["per_pixel"].groupby("dot_name"):
		axes2.plot(subframe.loc[dn,"mean"],color="grey",alpha=0.5,lw=1)
	#axes2.fill_between(mean.index,mean-std,mean+std,color="grey",alpha=0.2)
	axes2.plot(mean,color="k",lw=3)
	axes2.set(ylabel="Brightness")
	axes2.set_xticks(np.arange(1,13,3))
	axes2.set_xticklabels(["Jan","Apr","Jul","Oct"])
	axes2.set_yticks([])

	## Save it
	plt.savefig("..\\_plots\\transmission_vs_lights_inset.pdf")

	plt.show()