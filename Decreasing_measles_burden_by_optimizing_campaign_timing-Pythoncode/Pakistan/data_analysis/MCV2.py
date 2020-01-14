""" Does Pakistan have significant mcv2?"""
import sys
sys.path.insert(0, "..\\..\\")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import shapefile
from riskmap3.map_maker import *
from riskmap3.data_process.shapefile_tools import NamesFromSf
from riskmap3.data_process.linelist_tools import *
from riskmap3.data_process.timeseries_tools import *
import pickle

def CreateDotNames(df):
	df["dot_name"] = "asia:pakistan:"+df.province.str.lower()
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
	rejected = df[df.classification == 0][["time","dot_name","measles igm results","num_vaccine_doses"]]

	## Clean up the num_vaccine doses catagory
	clean_up = {"two":2,"Two":2,"zero":0,"Zero":0,"one":1,9:np.nan,"DK":np.nan,"Unknown":np.nan,"Nil":0,
				"One":1,"ONE":1,"UK":np.nan,"Three":3,"TWO":2,"NIL":0,"five":5,"Four":4}
	rejected["num_vaccine_doses"] = rejected.num_vaccine_doses.apply(lambda x: clean_up.get(x,x))
	
	## Drop NANs
	rejected = rejected[rejected.num_vaccine_doses.notnull()]

	## Create by state estimates
	mcv = rejected.groupby("dot_name")["num_vaccine_doses"].value_counts()
	total_entries = rejected.groupby("dot_name")["num_vaccine_doses"].count()
	mcv = (mcv/total_entries).sort_index(level=["dot_name","num_vaccine_doses"]).rename("mcv coverage")
	
	## MCV1 for comparison with DHS
	#mcv1 = mcv.loc["asia:pakistan:punjab"]
	mcv1 = mcv.groupby(level=0).apply(lambda x: np.sum(x.as_matrix()[1:])).rename("mcv1")
	mcv2 = mcv.groupby(level=0).apply(lambda x: np.sum(x.as_matrix()[2:])).rename("mcv2")

	## Make a Space time df for mcv2
	time = pd.DatetimeIndex(start="12/31/2008",end="12/31/2020",freq="SM")
	ST_df = EmptySpaceTimeDf(dot_names,time).rename("mcv2")
	for dot_name in dot_names:
		ST_df.loc[dot_name,:] = mcv2.loc[dot_name]
	
	## Pickle the results
	ST_df.to_pickle("..\\pickle_jar\\extrapolated_mcv2.pkl")

	## For comparison
	dhs = pd.read_pickle("..\\pickle_jar\\extrapolated_ri.pkl").rename("DHS coverage")
	dhs = dhs.groupby("admin").mean().sort_index()
	print(mcv1)
	print(dhs)
	print(mcv2)
