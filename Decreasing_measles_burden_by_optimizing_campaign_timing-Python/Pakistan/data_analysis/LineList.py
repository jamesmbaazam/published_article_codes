""" Functions to process the merged linelist data at admin1."""

## Append the location of the riskmap3 lib
import sys
sys.path.insert(0, "..\\..\\")

## Standard imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## For plotting, aliasing, etc.
import shapefile
from riskmap3.map_maker import *
from riskmap3.data_process.shapefile_tools import NamesFromSf
from riskmap3.data_process.linelist_tools import *
plt.rcParams["font.size"] = 24

## For saving the alias table
import pickle



def CreateDotNames(df):
	df["dot_name"] = "asia:pakistan:"+df.province.str.lower()
	return df

def Bin(df):
	t_index = pd.DatetimeIndex(start="12/31/2008",end="12/31/2017",freq="SM")
	data = {}
	for admin, sf in df.groupby("dot_name"):
		labeled = pd.cut(sf.time,t_index,labels=t_index[1:])
		data[admin] = labeled.value_counts().reindex(t_index[1:]).fillna(0.)

	## Create a data frame with the appropriate names	
	new_df = pd.concat(data.values(), keys=data.keys())
	new_df.index.rename(["admin","time"],inplace=True)
	return new_df

if __name__ == "__main__":

	## Where does the data live?
	data_dir = "..\\_data\\CaseData\\"

	## Get the merged dataset
	df = pd.read_pickle(data_dir+"merged_dataset.pkl")

	## And the shapefile
	## Get the shapefile and create a sf object
	shp = "..\\_data\\Shapefile\\Pakistan.shp"
	sf = shapefile.Reader(shp)
	dot_names = NamesFromSf(sf,admin_level=1)

	## Make the dot name column and keep records only
	## if they have a dot_name
	df = CreateDotNames(df)
	df = df[df.dot_name.notnull()]

	## Standardize the dot_name column
	## This is done with a pickled alias table, but setting table to initially be
	## empty (i.e. table = {}) will let you try the interactive alias experience.
	interactive_alias = False
	if interactive_alias:
		table = {}
	else:
		table = pickle.load(open(data_dir+"admin1_match_table.pkl","rb"))
	df, table, unused = InteractiveStandardizeNames(dot_names,df,tol=97.,match_table=table,
												drop_col=True)
	with open(data_dir+'admin1_match_table.pkl', 'wb') as file:
			pickle.dump(table, file, protocol=pickle.HIGHEST_PROTOCOL)

	## Create a time column
	df["time"] = df.date_onset
	df = df[df.time.notnull()]

	## Parse the different measels igm catagories
	catagories = {"POSITIVE":1,"Measles Lab-confirmed":1,"Measles Epi-Linked":1}
	df["classification"] = df["measles igm results"].apply(lambda s: catagories.get(s,0))

	## Slice to find cases and rejected cases
	cases = df[df.classification == 1][["time","dot_name","measles igm results"]]
	rejected = df[df.classification == 0][["time","dot_name","measles igm results"]]

	## And print some diagnostics by year...
	print("Total entries by year...")
	print(df.groupby(df.time.apply(lambda t: t.year))["dot_name"].count())
	print("Total confirmed by year...")
	print(cases.groupby(cases.time.apply(lambda t: t.year))["dot_name"].count())
	print("Total rejected by year...")
	print(rejected.groupby(rejected.time.apply(lambda t: t.year))["dot_name"].count())

	## Create the monthly time index you want
	cases = Bin(cases).rename("cases")
	rejected = Bin(rejected).rename("rejected")
	
	## Smooth the case data
	cases_smooth = cases.groupby("admin").rolling(window=2).mean().reset_index(level=0,drop=True)

	## Make pickles
	cases.to_pickle("..\\pickle_jar\\cases.pkl")
	rejected.to_pickle("..\\pickle_jar\\rejected.pkl")
	cases_smooth.to_pickle("..\\pickle_jar\\smoothed_cases.pkl")
	
	## Make a plot
	province = "asia:pakistan:punjab"
	fig, axes = plt.subplots(2,1,sharex=True,figsize=(12,8))
	axes[0].plot(cases.loc[province],label="Cases")
	axes[0].plot(rejected.loc[province],label="Rejected")
	axes[0].plot(cases_smooth.loc[province],label="Smoothed")
	axes[0].legend()
	axes[0].set(title=province[len("asia:pakistan:"):].title())
	
	## National level
	n_cases = cases.groupby("time").sum()
	n_rejected = rejected.groupby("time").sum()
	n_cases_smooth = cases_smooth.groupby("time").sum()
	axes[1].plot(n_cases)
	axes[1].plot(n_rejected)
	axes[1].plot(n_cases_smooth)
	axes[1].set(title="National level")

	## Aggregated over years
	p_name = province[len("asia:pakistan:")].upper()+province[len("asia:pakistan:")+1:]
	by_month = n_cases_smooth.groupby(lambda t: t.month).sum()
	by_month_p = cases_smooth.loc[province].groupby(lambda t: t.month).sum()
	fig, axes = plt.subplots(figsize=(12,6))
	axes.bar(by_month.index,by_month,width=0.75,color="grey",alpha=0.75,label="National")
	axes.bar(by_month_p.index,by_month_p,width=0.75,color="C0",label=p_name)
	axes.legend()
	axes.set(xlabel="Month",ylabel="Cases (confirmed)",title="Pakistan Seasonality")


	rejected_smoothed = rejected.groupby("admin").rolling(window=2).mean().reset_index(level=0,drop=True)
	n_rejected_smooth = rejected_smoothed.groupby("time").sum().loc["01-01-2010":]
	fig, axes = plt.subplots(figsize=(16,6))
	axes.grid(color="grey",alpha=0.2)
	axes.fill_between(n_rejected_smooth.index,len(n_rejected_smooth)*[0],n_rejected_smooth.values,color="grey",alpha=0.3)
	axes.plot(n_rejected_smooth,color="k",lw=2)
	axes.set(ylabel="Lab rejected cases",ylim=(0.,None))
	plt.tight_layout()
	plt.savefig("..\\_plots\\rejected_timeseries.pdf")
	plt.show()