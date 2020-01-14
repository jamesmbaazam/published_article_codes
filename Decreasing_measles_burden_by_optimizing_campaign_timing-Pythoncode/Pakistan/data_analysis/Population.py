import sys
sys.path.insert(0, "..\\..\\")

from datetime import datetime
starttime = datetime.now()

import numpy as np
import pandas as pd

## For shape files
import shapefile

## For tiff files
import georaster
import matplotlib.pyplot as plt
import matplotlib.patches as patches

## Risk map tools
from riskmap3.map_maker import *
from riskmap3.data_process.shapefile_tools import *
from riskmap3.data_process.timeseries_tools import *
plt.rcParams["font.size"] = 24.

## For interpolation
from scipy.interpolate import griddata

def WorldPopSeries(admin_level,sf,pop_file,data_dir,check_pickles=True,save=True,interpolate=[]):

	"""Function to create a series of pop estimates at admin level based on the
	world pop estimate. 

	check_pickles: bool, telling the function to check the data dir for premade versions
				   to save time (since computing them takes a while).

	save: save a pickle if you compute something new."""

	## Based on the naming convention, this is the pickle to
	## check (and save to)
	pickle = data_dir+"WorldPop\\_processed\\admin"+str(admin_level)+"_pop_"+pop_file[-9:-4]+".pkl"
	if check_pickles:
		try:
			population = pd.read_pickle(pickle)
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
		image = georaster.SingleBandRaster(data_dir+pop_file)

		## Remove zero pixels
		mask = image.r > 0.

		## Get GPS and pop, then delete the 
		## georaster to free up memory for the 
		## data frame.
		print("Removing empty pixels...")
		X, Y = image.coordinates()
		pop = image.r[mask]
		X = X[mask]
		Y = Y[mask]
		del image	

		## Create a data frame with admin level column.
		print("Creating a data frame...")
		df = np.array([np.arange(len(X),dtype=int),X,Y,pop]).T
		df = pd.DataFrame(df, columns=["pt","X","Y","pop"])
		del X, Y, pop
		print("Assigning admin"+str(admin_level)+" labels...")
		df = AdminLevelFromGPS(df,sf,admin_level,col_cluster="pt",col_x="X",col_y="Y",output_col="admin")

		## Compute the population, and release the 
		## data frame.
		print("Computing total population...")
		population = df.groupby("admin")["pop"].sum()
		del df

		## Save this run (before interpolation, so you can experiment with
		## different interpolation methods if needed).
		if save:
			print("Saving as "+pickle)
			population.to_pickle(pickle)
	
	## Interpolate missing parts of the tiff if needed, done regardless
	## of the pickle incase you want to overwrite nan's in different ways
	## and compare.
	names = NamesFromSf(sf,admin_level)
	if interpolate and (len(population.index) != len(names)):
		
		print("Interpolating missing admin units...")
		
		## Get centers
		admin_units = AdminLevelSubset(sf, admin_level)
		centers = {n:ShapeCenter(s) for n,s in admin_units.items()}
		centers = pd.DataFrame(list(centers.values()), 
							   index=pd.Index(centers.keys(),name="admin_level"),
							   columns=["center_x","center_y"])

		## Merge the data frames to force nan's at 
		## missing admin units
		population = pd.concat([centers, population],axis=1)
		
		## Now interpolate
		points = data[["center_x","center_y"]]
		z = data[["pop"]]
		for method in interpolate:
			population[method] = griddata(points,z,population[["center_x","center_y"]],method=method)
			population["pop"].fillna(population[method],inplace=True)
		population = population["pop"]

	return population

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

	## Set the admin level. We'll do it at 
	## admin 2 since that's more general and then upsample
	## afterwards.
	admin_level = 2

	## Get the shapefile and create a sf object
	shp = "..\\_data\\Shapefile\\Pakistan.shp"
	sf = shapefile.Reader(shp)

	## Get the world pop tiff
	data_dir = "..\\_data\\"
	pop_files = ["Worldpop\\popmap10adj.tif","Worldpop\\popmap15adj.tif"]
	years = [pd.to_datetime(2010,format="%Y"),pd.to_datetime(2015,format="%Y")]

	## Get the series
	series = []
	for file in pop_files:
		start_time = datetime.now()
		print("\nStarting computation for: "+file)
		pop = WorldPopSeries(admin_level,sf,file,data_dir,interpolate=["linear"])
		print("... done in {}".format(datetime.now()-start_time))
		series.append(pop)
	
	## Create multi index df
	df = SpaceTimeDf(series,years)

	## Interpolate over time
	freq = "SM"
	time = pd.DatetimeIndex(start="12/31/2008",end="12/31/2020",freq=freq)
	df = TimeInterpolate(df,freq,time)

	## Show results
	plot_time = "2015-11-30"
	population = df[:,plot_time]

	## Fig 1: a population map
	fig, axes = plt.subplots(figsize=(12,8))
	PlotDfonSf(fig,axes,population,sf,admin_level={2},alpha=1.,colorbar=True,
			   cmap=matplotlib.cm.viridis,missing_color="red",vlim=(0,4000000))
	PlotBorders(fig,axes,sf,admin_level={0})
	axes.set(title="Population, "+plot_time[:plot_time.find("-")],aspect="equal")
	axes.axis("off")
	plt.tight_layout()
	plt.savefig("..\\_plots\\pop_map.png")

	## Make an admin2 pickle
	df = PortToShapefileEdit(df)
	#df.to_pickle("..\\pickle_jar\\population_admin2.pkl")

	## Upsample to admin 1
	#df = df.groupby([lambda s: s[0][:s[0].rfind(":")],"time"]).sum()

	## Pickle the result
	#df.to_pickle("..\\pickle_jar\\extrapolated_population.pkl")
	
	## Plot a time series
	admin = "asia:pakistan:punjab:lahore"
	admin_name = "Lahore, Punjab"
	timeseries = df.loc[admin]
	raster_values = [s.loc[admin] for s in series] 
	fig, axes = plt.subplots(figsize=(10,6))
	axes.grid(color="grey",alpha=0.2)
	axes.plot(timeseries,color="k",label="Population estimate")
	axes.plot(years,raster_values,marker="o",ls="None",markersize=12,
			  color="xkcd:red wine",label="Raster values")
	axes.legend()
	axes.set(title=admin_name)
	axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.tight_layout()
	plt.savefig("..\\_plots\\pop_series_example.png")
	plt.show()
