""" UpdatedLineList.py

On 7/5/2018, we recieved an updated 2017 line list document from Quamrul.
This script is processing for that particular document (and 2018 when it comes) to
be compared to model projections from the original dataset. """
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
plt.rcParams["font.size"] = 28.


## For saving the alias table
import pickle



def CreateDotNames(df):
	df["dot_name"] = "asia:pakistan:"+df.Province.str.lower()+":"+df.Districts.str.lower()
	return df

def ParseDate(df,col_name,errors="coerce"):

	""" This is a wildly inefficient way to do this, but need something that works now."""

	## First correct those in non-standard but interpretable formats
	def preprocess(t):
		if isinstance(t,pd.datetime):
			return t
		t = str(t)
		t = t.replace(".","").replace("Z","")
		if t.find("/") == -1:
			return t[:2]+"/"+t[2:4]+"/"+t[4:]
		else:
			return t
	df[col_name] = [pd.to_datetime(preprocess(t),errors=errors) for t in df[col_name]]
	return df

def Bin(df):
	t_index = pd.DatetimeIndex(start="12/31/2016",end="12/31/2017",freq="SM")
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

	## Get the file
	df = pd.read_excel(data_dir+"_raw\\Measles Januray-December 2017.xlsx",
					   sheetname="Sheet1",header=3,
					   na_values=["-"])

	## And the shapefile
	## Get the shapefile and create a sf object
	shp = "..\\_data\\Shapefile\\Pakistan_edit.shp"
	sf = shapefile.Reader(shp)
	dot_names = NamesFromSf(sf,admin_level=2)

	## Make the dot name column and keep records only
	## if they have a dot_name
	df = CreateDotNames(df)
	df["dot_name"] = df.dot_name.str.replace("pakistan:kpk","pakistan:khyber pakhtoon")
	df["dot_name"] = df.dot_name.str.replace("pakistan:kp","pakistan:khyber pakhtoon")
	df["dot_name"] = df.dot_name.str.replace("pakistan:g.b","pakistan:gilgit baltistan")
	df["dot_name"] = df.dot_name.str.replace("pakistan:gb","pakistan:gilgit baltistan")
	df = df[df.dot_name.notnull()]
	
	## Standardize the dot_name column
	## This is done with a pickled alias table, but setting table to initially be
	## empty (i.e. table = {}) will let you try the interactive alias experience.
	interactive_alias = False
	if interactive_alias:
		table = {}
	else:
		table = pickle.load(open(data_dir+"admin2_match_table.pkl","rb"))
	df, table, unused = InteractiveStandardizeNames(dot_names,df,tol=97.,match_table=table,
												drop_col=True)
	with open(data_dir+'admin2_match_table.pkl', 'wb') as file:
			pickle.dump(table, file, protocol=pickle.HIGHEST_PROTOCOL)

	## Parse the time column
	df = ParseDate(df,"Date of Onset")
	df["time"] = df["Date of Onset"]
	df = df[df.time.notnull()]

	## Parse the different measels igm catagories
	catagories = {"POSITIVE":1}
	df["classification"] = df["Measles IgM Results"].apply(lambda s: catagories.get(s,0))

	## Pickle the minimally processed line list
	df.to_pickle("..\\_data\\CaseData\\updated2017.pkl")

	## Slice to find cases and rejected cases
	admin2_cases = df[df.classification == 1][["time","dot_name"]]

	## Create the monthly time index you want
	cases = Bin(admin2_cases).rename("cases")

	## Upsample to Admin1
	cases = cases.groupby([lambda s: s[0][:s[0].rfind(":")], "time"]).sum()
	cases.index.rename(["admin","time"],inplace=True)

	## Subset to the right provinces (i.e. not FATA and AJK)
	provinces = ["asia:pakistan:"+p for p in ["punjab","sindh","balochistan","islamabad","khyber pakhtoon","gilgit baltistan"]]
	cases = cases.loc[provinces,:]
	
	## Smooth the case data
	cases_smooth = cases.groupby("admin").rolling(window=2).mean().reset_index(level=0,drop=True).dropna()
	
	## Save a pickle
	cases.to_pickle("..\\pickle_jar\\updated_cases.pkl")
	cases_smooth.to_pickle("..\\pickle_jar\\updated_smoothed_cases.pkl")
	print(cases_smooth)

	## Get the old dataset for comparison
	og_cases = pd.read_pickle("..\\pickle_jar\\cases.pkl").rename("og_cases")
	og_smooth = pd.read_pickle("..\\pickle_jar\\smoothed_cases.pkl").rename("og_smooth")
	
	## National level comparison
	national = cases.groupby("time").sum()
	nat_smooth = cases_smooth.groupby("time").sum()
	og_national = og_cases.groupby("time").sum()
	og_smooth = og_smooth.groupby("time").sum()

	## Plot it up
	fig, axes = plt.subplots(figsize=(12,8))
	axes.grid(color="grey",alpha=0.2)
	axes.plot(og_national, color="k", alpha=0.7, lw=3,label="Original dataset")
	axes.plot(national,color="C3", label="Updated dataset")
	axes.legend()
	axes.set(ylabel="Reported cases")
	plt.tight_layout()
	
	fig, axes = plt.subplots(figsize=(12,8))
	axes.grid(color="grey",alpha=0.2)
	axes.plot(og_smooth, color="k", alpha=0.7, lw=3,label="Original dataset")
	axes.plot(nat_smooth,color="C0", label="Updated dataset")
	axes.legend()
	axes.set(ylabel="Smoothed cases")
	plt.tight_layout()

	## What about cases per 100k?
	total_cases = cases.groupby(level=0).sum()
	population = pd.read_pickle("..\\pickle_jar\\extrapolated_population.pkl")
	population = population.loc[:,"01-01-2017":"12-31-2017"].groupby(level=0).mean()
	case_density = (1.e5*total_cases/population).rename("Cases per 100k in 2017")
	print(case_density)

	### What's going on in KP?
	kp = admin2_cases.loc[admin2_cases.dot_name.str.contains("khyber pakhtoon")]

	## Split into pre and post September 2017
	cut_date = pd.to_datetime("08-01-2017",format="%m-%d-%Y")
	mask = (kp.time <= cut_date) & (kp.time > pd.to_datetime("05-01-2017",format="%m-%d-%Y"))
	kp_pre = kp[mask]["dot_name"].value_counts()
	kp_post = kp[kp.time > cut_date]["dot_name"].value_counts()

	## Make some maps
	province_sf = [sr.shape for sr in sf.iterShapeRecords() if sr.record[1].lower() == "asia:pakistan:khyber pakhtoon"][0]
	bbox = province_sf.bbox
	fig, axes = plt.subplots(figsize=(10,8))
	PlotDfonSf(fig,axes,kp_pre,sf,admin_level={2},alpha=1.,colorbar=True,
			   cmap=matplotlib.cm.Oranges,missing_color="0.85")
	PlotBorders(fig,axes,sf,admin_level={0,2})
	axes.set(aspect="equal",
			 xlim=(0.99*bbox[0],1.01*bbox[2]),
			 ylim=(0.99*bbox[1],1.01*bbox[3]),
			 title="Pre "+str(cut_date))
	axes.axis("off")
	plt.tight_layout()

	fig, axes = plt.subplots(figsize=(10,8))
	PlotDfonSf(fig,axes,kp_post,sf,admin_level={2},alpha=1.,colorbar=True,
			   cmap=matplotlib.cm.Oranges,missing_color="0.85")
	PlotBorders(fig,axes,sf,admin_level={0,2})
	axes.set(aspect="equal",
			 xlim=(0.99*bbox[0],1.01*bbox[2]),
			 ylim=(0.99*bbox[1],1.01*bbox[3]),
			 title="Post "+str(cut_date))
	axes.axis("off")
	plt.tight_layout()


	## Why do we need to smooth?
	combined = pd.concat([og_national.loc[:"12-31-2016"],national],axis=0).loc["01-01-2012":]
	combined_smooth = pd.concat([og_smooth.loc[:"12-31-2016"],nat_smooth],axis=0).loc["01-01-2012":]
	fig, axes = plt.subplots(figsize=(16,6))
	#axes.fill_between(combined.index,len(combined)*[0.],combined.values,
	#				  color="#FF420E",alpha=0.2)
	axes.fill_between(combined_smooth.index,len(combined_smooth)*[0.],combined_smooth.values,
					  color="#68829E",alpha=0.3)
	axes.plot(combined_smooth,lw=2,color="#68829E",label="Smoothed case data")
	axes.plot(combined,color="#FF420E",label="Case data")
	handles, labels = axes.get_legend_handles_labels()
	axes.legend(handles[::-1],labels[::-1],loc=1)
	axes.set(ylim=(0,None),title="Lab confirmed measles in Pakistan")
	plt.tight_layout()
	plt.savefig("..\\_plots\\case_data_overview.png")





	plt.show()
