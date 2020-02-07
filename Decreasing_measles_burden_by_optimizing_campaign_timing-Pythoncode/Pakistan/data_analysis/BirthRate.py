""" BirthRate.py
Birth rate estimation from the Worldpop birth data."""
from datetime import datetime
starttime = datetime.now()

## Set up path to risk map tools
import sys
sys.path.insert(0,"..\\..\\")

import numpy as np
import pandas as pd

## For shape files
import shapefile

## For tiff files
import georaster
import matplotlib.pyplot as plt

## Risk map tools
from riskmap3.map_maker import *
from riskmap3.data_process.shapefile_tools import *
from riskmap3.data_process.timeseries_tools import *

## For interpolation
from scipy.interpolate import griddata

def WorldPopSeries(admin_level,sf,file,data_dir,check_pickles=True,save=True,interpolate=[]):

	"""Function to create a series of pop estimates at admin level based on the
	world pop estimate. 

	check_pickles: bool, telling the function to check the data dir for premade versions
				   to save time (since computing them takes a while).

	save: save a pickle if you compute something new."""

	## Based on the naming convention, this is the pickle to
	## check (and save to)
	i = file.find("PAK")
	pickle = data_dir+"BirthsWorldPop\\_processed\\admin"+str(admin_level)+"_births_"+file[i:i+7]+".pkl"
	if check_pickles:
		try:
			births = pd.read_pickle(pickle)
			compute = False
		except:
			compute = True
	else:
		compute = True

	## If we get here, we have to compute it.
	if compute:
		print("No pickle matching this request ("+pickle+")! Starting computation...")

		## Get the tif as a georaster
		print("Getting the georaster...")
		image = georaster.SingleBandRaster(data_dir+file)

		## Remove zero pixels
		mask = image.r > 0.

		## Get GPS and births, then delete the 
		## georaster to free up memory for the 
		## data frame.
		print("Removing empty pixels...")
		X, Y = image.coordinates()
		births = image.r[mask]
		X = X[mask]
		Y = Y[mask]
		del image, mask
		
		## Create a data frame with admin level column.
		print("Creating a data frame...")
		df = np.array([np.arange(len(X),dtype=int),X,Y,births]).T
		df = pd.DataFrame(df, columns=["pt","X","Y","births"])
		del X, Y, births
		print("Assigning admin"+str(admin_level)+" labels...")
		df = AdminLevelFromGPS(df,sf,admin_level,col_cluster="pt",col_x="X",col_y="Y",output_col="admin")

		## Compute the births, and release the 
		## data frame.
		print("Computing total population...")
		births = df.groupby("admin")["births"].sum()
		del df

		## Save this run (before interpolation, so you can experiment with
		## different interpolation methods if needed).
		if save:
			print("Saving as "+pickle)
			births.to_pickle(pickle)
	
	## Interpolate missing parts of the tiff if needed, done regardless
	## of the pickle incase you want to overwrite nan's in different ways
	## and compare.
	names = NamesFromSf(sf,admin_level)
	if interpolate and (len(births.index) != len(names)):
		
		print("Interpolating missing admin units...")
		
		## Get centers
		admin_units = AdminLevelSubset(sf, admin_level)
		centers = {n:ShapeCenter(s) for n,s in admin_units.items()}
		centers = pd.DataFrame(centers.values(), 
							   index=pd.Index(centers.keys(),name="admin_level"),
							   columns=["center_x","center_y"])

		## Merge the data frames to force nan's at 
		## missing admin units
		births = pd.concat([centers, births],axis=1)
		data = births.dropna()
		
		## Now interpolate
		points = data[["center_x","center_y"]]
		z = data[["births"]]
		for method in interpolate:
			births[method] = griddata(points,z,births[["center_x","center_y"]],method=method)
			births["births"].fillna(births[method],inplace=True)
		births = births["births"]

	return births

def SpatialInterpolate(admin_level,sf,series,methods=["linear"]):

	## Interp variable
	var = series.name
		
	## Get centers
	admin_units = AdminLevelSubset(sf, admin_level)
	centers = {n:ShapeCenter(s) for n,s in admin_units.items()}
	centers = pd.DataFrame(centers.values(), 
							   index=pd.Index(centers.keys(),name="admin_level"),
							   columns=["center_x","center_y"])

	## Merge the data frames to force nan's at 
	## missing admin units
	series = pd.concat([centers, series],axis=1)
	data = series.dropna()
		
	## Now interpolate
	points = data[["center_x","center_y"]]
	z = data[[var]]
	for method in methods:
		series[method] = griddata(points,z,series[["center_x","center_y"]],method=method)
		series[var].fillna(series[method],inplace=True)
	series = series[var]

	return series

def PortToShapefileEdit(df):

	""" For admin-2, the Karachi portion of the shapefile needed to be edited. This function
	implements that edit."""

	## If we want the df to align with the edited shape file, we have to correct
	## for karachi.
	karachi = df.loc[df.index.get_level_values(level=0).str.startswith("asia:pakistan:sindh:khi")]
	karachi = karachi.groupby(level=1).sum()
	karachi = pd.Series(karachi.values,index=[["asia:pakistan:sindh:karachi"]*len(karachi),karachi.index])
	karachi.index.rename(["admin","time"],inplace=True)

	## Now drop the old karachi districts (the split ones)
	## and add the new karachi district.
	df = df.loc[~df.index.get_level_values(level=0).str.startswith("asia:pakistan:sindh:khi")]
	df = pd.concat([df,karachi],axis=0)

	return df

if __name__ == "__main__":

	## Get the shapefile and create a sf object
	admin_level = 2
	shp = "..\\_data\\Shapefile\\Pakistan.shp"
	sf = shapefile.Reader(shp)

	## Choose files to use
	data_dir = "..\\_data\\"
	files =	["BirthsWorldPop\\PAK2010adjustedBirths.tif",
			 "BirthsWorldPop\\PAK2012adjustedBirths.tif",
			 "BirthsWorldPop\\PAK2015adjustedBirths.tif",
			 "BirthsWorldPop\\PAK2020adjustedBirths.tif"]
	years = [pd.to_datetime(t) for t in ["2010","2012","2015","2020"]]

	## Process tiffs
	series = []
	for file in files:
		births = WorldPopSeries(admin_level,sf,file,data_dir)
		print(datetime.now() - starttime)
		startime = datetime.now()
		series.append(births)

	## Make a space time df
	time = pd.DatetimeIndex(start="12/31/2008",end="12/31/2020",freq="SM")
	births = SpaceTimeDf(series,years)
	births = TimeInterpolate(births,freq="SM",time=time)

	## Make an admin2 births
	births = PortToShapefileEdit(births)
	births.to_pickle("..\\pickle_jar\\births_admin2.pkl")

	## Upsample it
	births = births.groupby([lambda s: s[0][:s[0].rfind(":")], "time"]).sum()

	## Get the population df
	population = pd.read_pickle("..\\pickle_jar\\extrapolated_population.pkl")
	
	## Convert to birthrate
	birth_rate = 1000.*births/population
	birth_rate = 1000.*((1. + birth_rate/1000.)**(1./24.) - 1.)
	
	## Pickle the result
	birth_rate.to_pickle("..\\pickle_jar\\extrapolated_birth_rate.pkl")
	
	## plots
	## use the dhs for comparison...
	dhs = pd.read_pickle("..\\pickle_jar\\extrapolated_dhs_birth_rate.pkl")
	dhs_avg = dhs.groupby(level=0).mean().rename("DHS")
	births_avg = births.groupby(level=0).mean().rename("Births")
	pop_avg = population.groupby(level=0).mean().rename("Population")
	
	## Map
	fig, axes = plt.subplots(figsize=(12,10))
	avg = birth_rate.groupby(level=0).mean().rename("WorldPop")
	print(pd.concat([avg,dhs_avg,pop_avg,births_avg],axis=1))
	PlotDfonSf(fig,axes,avg,sf,admin_level={1},alpha=1.,colorbar=True,
			   cmap=matplotlib.cm.viridis)
	PlotBorders(fig,axes,sf,admin_level={0,1})
	axes.set(aspect="equal",title="Birthrate (WorldPop)")
	axes.axis("off")
	plt.tight_layout()

	## Time series
	province = "asia:pakistan:sindh"
	fig, axes = plt.subplots(figsize=(12,6))
	axes.plot(birth_rate.loc[province],label="WorldPop",c="C3")
	try:
		axes.plot(dhs.loc[province],label="DHS",c="k")
	except:
		axes.plot([],label="DHS (no clusters)",c="k")
	axes.legend()
	axes.set(title="Birth rate over time, "+province[province.rfind(":")+1:].title())
	plt.show()
	