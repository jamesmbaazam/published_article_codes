""" Exploration of infection age in the line list data."""
from __future__ import print_function
import sys
sys.path.insert(0, "..\\..\\")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import shapefile
from riskmap3.map_maker import *
from riskmap3.data_process.shapefile_tools import NamesFromSf
from riskmap3.data_process.linelist_tools import *
import pickle
plt.rcParams["font.size"] = 28.

def CreateDotNames(df):
	df["dot_name"] = "asia:pakistan:"+df.province.str.lower()#+":"+df.district.str.lower()
	return df

if __name__ == "__main__":

	data_dir = "..\\_data\\CaseData\\"

	## Get the merged dataset
	df = pd.read_pickle(data_dir+"merged_dataset.pkl")

	## And the shapefile
	## Get the shapefile and create a sf object
	shp = "..\\_data\\Shapefile\\Pakistan.shp"
	sf = shapefile.Reader(shp)
	dot_names = NamesFromSf(sf,admin_level=1)

	## Make the dot name column
	df = CreateDotNames(df)
	df = df[df.dot_name.notnull()]

	## Standardize the dot_name column
	table = pickle.load(open(data_dir+"admin1_match_table.pkl","rb"))
	df, table, unused = InteractiveStandardizeNames(dot_names,df,tol=97.,match_table=table,
												drop_col=True)

	## Create a time column
	df["time"] = df.date_onset
	df = df[df.time.notnull()]

	## Parse the different measels igm catagories
	catagories = {"POSITIVE":1,"Measles Lab-confirmed":1,"Measles Epi-Linked":1}
	df["classification"] = df["measles igm results"].apply(lambda s: catagories.get(s,0))

	## Slice to find cases and rejected cases
	cases = df[df.classification == 1][["time","dot_name","age_months"]]

	## Clean up
	cases = cases[cases.age_months.notnull()]
	cases = cases[cases.age_months <= 720.] ## No one over 60 years old.

	## Create histogram bins
	bins = np.arange(0.,721.,12)
	width = 12.

	## Groupby, cut, and reshape
	histograms = {}
	for state, subframe in cases.groupby("dot_name"):
		binned = pd.cut(subframe.age_months,bins,labels=bins[1:])
		hist = binned.value_counts().sort_index().rename("infection_age")
		histograms[state] = hist
	inf_age = pd.concat(histograms.values(),keys=histograms.keys())
	inf_age.index.rename(["admin","age"],inplace=True)

	## Create a national level histogram
	national = inf_age.groupby("age").sum()
	
	## Select a provinces for analysis too
	provinces = ["asia:pakistan:punjab","asia:pakistan:sindh"]
	p_names = [province[province.rfind(":")+1:].title() for province in provinces]
	dfs = [inf_age.loc[province] for province in provinces]
	colors = ["C"+str(i % 9) for i in range(len(provinces))]

	## Accumulate the dfs (so the plot is stacked)
	adfs = []
	base = 0.*dfs[0]
	for df in dfs:
		base = base + df
		adfs.append(base)

	## Make some plots
	fig, axes = plt.subplots(figsize=(12,8))
	axes.bar(national.index,national,width=width,color="grey",label="National")
	for c, p_name, df in zip(colors,p_names[::-1],adfs[::-1]):
		axes.bar(df.index,df,width=width,color=c,label=p_name)
	axes.legend()
	axes.set(title="Infection ages",xlabel="Age (months)")#,xlim=(0.,180.))
	

	## A simpler plot
	bins = range(0,181,6)
	national = pd.cut(cases.age_months,bins).value_counts().sort_index().rename("count")
	national = pd.concat([national,pd.Series(range(len(national)),index=national.index,name="bin_number")],axis=1)

	## Print some statistics
	total_cases = national["count"].sum()
	total_6m_to_5y = national.loc[pd.Interval(6,12):pd.Interval(54,60)]
	total_6m_to_10y = national.loc[pd.Interval(6,12):pd.Interval(114,120)]
	frac_6m_to_5y = total_6m_to_5y["count"].sum()/total_cases
	frac_6m_to_10y = total_6m_to_10y["count"].sum()/total_cases
	print("Fraction of cases in (6m, 5y] = {0:.3f}%".format(100.*frac_6m_to_5y))
	print("Fraction of cases in (6m, 10y] = {0:.3f}%".format(100.*frac_6m_to_10y))
	
	## Plot the full set
	w = 2./3.
	fig, axes = plt.subplots(figsize=(16,12))
	axes.axvspan(total_6m_to_10y["bin_number"].min()-0.75*w,total_6m_to_10y["bin_number"].max()+0.75*w,color="xkcd:red wine",alpha=0.4)
	axes.axvspan(total_6m_to_5y["bin_number"].min()-0.75*w,total_6m_to_5y["bin_number"].max()+0.75*w,color="xkcd:red wine",alpha=0.7)
	axes.bar(national["bin_number"].values,national["count"].values,width=w,color="k")

	## Add explanatory text
	axes.text(total_6m_to_5y["bin_number"].mean(),4150,"6 month to 5 year SIA\n"+r"{0:.0f}% of cases".format(100.*frac_6m_to_5y),
			  fontsize=28,horizontalalignment="center",verticalalignment="center",color="white")
	axes.text(1.45*total_6m_to_10y["bin_number"].mean(),4150,"6 month to 10 year SIA\n"+r"{0:.0f}% of cases".format(100.*frac_6m_to_10y),
			  fontsize=28,horizontalalignment="center",verticalalignment="center",color="white")

	## Set up the ticks
	axes.set_xticks(np.arange(0,len(national),2))
	axes.set_xticklabels(national.index[::2],rotation=45)
	axes.set_yticks([])
	axes.set_yticklabels([])

	## Remove some spines
	axes.spines["top"].set_visible(False)
	axes.spines["right"].set_visible(False)
	axes.spines["left"].set_visible(False)

	## Set up the axis labels, etc.
	axes.set(xlabel="Age bin (months)",ylim=(0,4500))
	axes.set_ylabel("Age distribution of cases",color="k")
	
	## Finish up
	plt.tight_layout()
	plt.savefig("..\\_plots\\age_distribution.pdf")

	plt.show()


