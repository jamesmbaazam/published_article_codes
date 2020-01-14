""" dhs_tools.py

Functions to facilitate working with DHS data.

As is, these functions are too Nigeria specific. We need to test these more generally."""

## Standard imports
import numpy as np
import pandas as pd

## For the GIS dbf location files
from dbfread import DBF

def GetData(path,data_dir="..\\..\\_data\\",columns=None,
			with_gps=False,geo_path=None,geo_columns=None):

	"""Use pandas built-in read stata function to extract DHS data located in
	data_dir. Path is the subfolder within data_dir and columns are the desired
	columns within the stata file (specifying subsets of the total columns will
	increase speed, sometimes dramatically).
	
	If with_gps, function uses geo_path to find a dbf file and merges it on cluster. If 
	geo_columns is not specficied, this is done with all columns in the dbf file."""

	## Check for compatability between options.
	if with_gps and geo_path is None:
		raise ValueError("Need a geo_path to the .dbf file to append spatial information.")

	## Read_stata for .dta files:
	if path.endswith(".DTA"):
		df = pd.read_stata(data_dir + path, columns=columns)

	## Otherwise throw an error.
	else:
		raise ValueError("Not a supported path type.")

	## Merge with gps data if required.
	if with_gps:
		geo_df = GetDBF(geo_path,data_dir)
		if geo_columns is None:
			geo_columns = geo_df.columns
		else:
			geo_columns = ["cluster"] + geo_columns
			geo_df = geo_df[geo_columns]

		## Merge with pd.merge
		df = df.merge(geo_df,left_on="v001",right_on="cluster")
		df.drop("cluster",axis=1,inplace=True)

	return df

def GetDBF(path,data_dir="..\\..\\_data\\"):

	"""Function to read a .dbf file and return a pandas dataframe."""

	if not path.endswith('.dbf'):
		raise ValueError("geo_path must be a .dbf file!")

	## Use dbfread to create a data frame.
	table = DBF(data_dir+path)
	df = pd.DataFrame(iter(table))
	
	## Add a cluster number field
	cluster_number = [int(s[6:]) for s in df["DHSID"]]
	df["cluster"] = cluster_number

	return df

if __name__ == "__main__":

	path = "DHS_RAW\\NIE_2008\\ngkr53dt\\NGKR53FL.DTA"
	geo_path = "DHS_RAW\\NIE_2008\\ngge52fl\\NGGE52FL.dbf"
	#path = "DHS_Raw\\NIE_2013\\ngkr6adt\\NGKR6AFL.DTA"
	#geo_path = "DHS_Raw\\NIE_2013\\ngge6afl\\NGGE6AFL.dbf"

	df = GetData(path,columns=["v001","hw1","h9"],with_gps=True,
				 geo_path=geo_path,geo_columns=["URBAN_RURA","LATNUM","LONGNUM"])
	print(df)

