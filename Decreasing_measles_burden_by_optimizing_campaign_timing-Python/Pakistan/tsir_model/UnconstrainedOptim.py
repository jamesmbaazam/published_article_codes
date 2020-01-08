""" UnconstrainedOptim.py

Optimizing the SIA timing without operational constraints, i.e. just trying one
every month of the year in 2018."""
import sys
sys.path.insert(0,"..\\..\\")

import numpy as np
import pandas as pd

## And then the standard
import matplotlib.pyplot as plt

## For plotting standardization
from riskmap3.map_maker import *

## Overwrite some risk map defaults
plt.rcParams["font.size"] = 30.

def describe_future(samples,I):
	future_cases_samples = samples[:,I:]
	total_samples = np.sum(future_cases_samples,axis=1)
	future_cases = total_samples.mean()
	#future_cases = np.percentile(total_samples,50.)
	future_cases_std = np.std(total_samples)
	low = np.percentile(total_samples, 2.5)
	high = np.percentile(total_samples, 97.5)
	return future_cases, future_cases_std, low, high

## Create the time index so I can slice the 
## samples appropriately.
time_index = pd.DatetimeIndex(start="15-01-2012",end="31-12-2020",freq="SM")

## Create an index cut-off corresponding to the appropriate time.
cut_off_time = pd.to_datetime("31-12-2017",format="%d-%m-%Y")
cut_off_index = np.argmin(np.abs(time_index - cut_off_time))

## Get the baseline (i.e. without SIA)
npz_file = np.load("_updated_npz2\\no_future.npz")
bl_mean, bl_std, _, _ = describe_future(npz_file["I_samples"],cut_off_index)

## Get the npz files and extract
## the infecteds sample summaries.
means = np.zeros((2,12))
stds = np.zeros((2,12))
lows = np.zeros((2,12))
highs = np.zeros((2,12))
direct_fraction = np.zeros((2,12))
susceptibles_immunized = np.zeros((2,12))

for k,y in enumerate(["_18","_19"]):
	for i in range(1,13):

		## Construct the file name
		fname = str(i)+y+"_0.npz"
		if len(fname[:fname.find("_")]) != 2:
			fname = "0"+fname

		## Get the SIA date
		sia_date = fname[:2]+"-20"+fname[3:5]+"-15"
		sia_date = pd.to_datetime(sia_date,format="%m-%Y-%d")

		## Pick a province?
		#fname = "islamabad_"+fname

		## Get the npz_file
		#npz_file = np.load("..\\npz_jar\\"+fname)
		npz_file = np.load("_updated_npz2\\"+fname)

		## Compute the summary stats
		mean, std, low, high = describe_future(npz_file["I_samples"],cut_off_index)

		## Estimate by month
		print("{} estimate = {} +/- {}".format(fname,mean,std))

		## Compute average S_t
		S_t = npz_file["S_samples"].mean(axis=0)

		## Find the index before and after the SIA
		i_before = np.argmin(np.abs(time_index - sia_date))
		i_after = np.argmin(np.abs(time_index - sia_date - pd.to_timedelta(30,unit="d")))

		## Compute the change in susceptibles
		change_in_S = S_t[i_before] - S_t[i_after]

		## Store the results
		means[k,i-1] = mean
		stds[k,i-1] = std
		lows[k,i-1] = low
		highs[k,i-1] = high
		direct_fraction[k,i-1] = change_in_S/(bl_mean-mean)
		susceptibles_immunized[k,i-1] = change_in_S

## Compute the averted infections
infections_averted = bl_mean - means

## Put the years together
means = means.reshape((24,))
stds = stds.reshape((24,))
lows = lows.reshape((24,))
highs = highs.reshape((24,))

## Find the minimum
i_min = np.argmin(means)

## Plot it up
m = np.arange(1,13)
fig, axes = plt.subplots(figsize=(8,8))
#axes.grid(color="grey",alpha=0.2)

## Plot the low seasons
axes.axvspan(5, 10, alpha=0.15, color="C0")
axes.axvspan(12+5, 12+10, alpha=0.15, color="C0")
#axes.plot([5,10],[4.55e6,4.55e6],lw=4,c="C0",alpha=0.5)
#axes.plot([12+5,12+10],[4.55e6,4.55e6],lw=4,c="C0",alpha=0.5)

## And the expected infections
#axes.plot(np.arange(1,13),means[:12],c="k",marker="s",markersize=10,ls="dashed")
#axes.plot(12+m,means[12:],c="C3",marker="o",markersize=10,ls="dashed")
#axes.errorbar(np.arange(1,13),means[:12],yerr=stds[:12],c="k",marker="s",markersize=10,ls="dashed")
#axes.errorbar(12+m,means[12:],yerr=stds[12:],c="C3",marker="o",markersize=10,ls="dashed")
#axes.plot(12+m[:5],means[12:17],c="C3",marker="o",markersize=10,ls="dashed")

## Fill between option
axes.fill_between(np.arange(1,13),means[:12]-stds[:12],means[:12]+stds[:12],color="grey",alpha=0.2)
axes.fill_between(12+m,means[12:]-stds[12:],means[12:]+stds[12:],color="C3",alpha=0.2)
#axes.fill_between(np.arange(1,13),lows[:12],highs[:12],color="grey",alpha=0.2)
#axes.fill_between(12+m,lows[12:],highs[12:],color="C3",alpha=0.2)
axes.plot(np.arange(1,13),means[:12],c="k",marker="s",markersize=10,ls="None")
axes.plot(12+m,means[12:],c="C3",marker="o",markersize=10,ls="None")

## Including the optimum
#axes.plot(np.arange(1,len(means)+1)[i_min],means[i_min],markeredgecolor="C4",marker="o",
#		  markersize=20,markerfacecolor="None",markeredgewidth=2)
axes.set(xlabel="SIA timing (month)",ylabel="Total infections (2018-21)")#,ylim=(None,2.4e6))
axes.ticklabel_format(axis="y",style="sci",scilimits=(0,1))
axes.set_xticks(np.arange(2,25,2))
axes.set_xlim((0,25))

## Some text for the low seasons
axes.text(7.5/25, 0.91, "2018 low\nseason",
		verticalalignment='center', horizontalalignment='center',
		color="k", fontsize=22,transform=axes.transAxes)
axes.text((12.+7.5)/25, 0.91, "2019 low\nseason",
		verticalalignment='center', horizontalalignment='center',
		color="C3", fontsize=22,transform=axes.transAxes)

axes.set_xticklabels(["2","4","6\n2018","8","10","12"]+["2","4","6\n2019","8","10","12"])
for xtick, color in zip(axes.get_xticklabels(), 6*["k"]+6*["C3"]):
	xtick.set_color(color)
#axes.legend()
#axes.set_ylim((None,5.e6))
plt.tight_layout()
#plt.savefig("..\\_plots\\unconstrained_updated.pdf")


## Make a stacked bar chart of direct fraction
herd_immunity = 1. - direct_fraction
fig, axes = plt.subplots(figsize=(16,10))
axes.bar(np.arange(1,13),direct_fraction[0,:],width=2./3.,color="#375E97",label="Vaccinated")
axes.bar(np.arange(1,13),herd_immunity[0,:],width=2./3.,bottom=direct_fraction[0,:],
		 color="#FB6542",label="Herd immunity")
axes.set_xticks(np.arange(1,13))
axes.set_xlim((0.5,12.5))
axes.set_xticklabels(["1","2","3","4","5","6","7","8","9","10","11","12"])
axes.set_yticks([])
axes.set_yticklabels([])
axes.set(xlabel="2018\nSIA timing (month)",title="Direct and indirect contibutions to averted infections\n",
		 ylabel="Fraction")
for i,h in enumerate(herd_immunity[0,:]):
	axes.text(i+1,1.,"{0:.2f}".format(h),
 			  horizontalalignment="center",verticalalignment="bottom",color="#FB6542")
axes.spines['right'].set_visible(False)
axes.spines['left'].set_visible(False)
axes.spines['top'].set_visible(False)
axes.legend(loc=4,framealpha=0.95)
plt.tight_layout()
plt.savefig("..\\_plots\\herd_immunity.pdf")



plt.show()


