""" RI_map.py

Script to use riskmap to make a map of routine immunization coverage. For now this is just from the summary
document, since there isn't available gps data for the DHS 2012-13"""

## Orient the script's path
## to include the riskmap3 lib
import sys
sys.path.insert(0,"..\\..\\")

## Standard stuff
import numpy as np
import pandas as pd

## For riskmapping
import shapefile
from riskmap3.map_maker import *
from riskmap3.data_process.shapefile_tools import *
from riskmap3.data_process.timeseries_tools import *
import matplotlib.patches as patches
plt.rcParams["font.size"] = 24.

def DHSSummary():

	""" Values for this study are taken directly from the DHS summary since the goal is
	admin level 1 results. """

	## The DHS years
	dhs_years = [pd.to_datetime("2006-06-15"),pd.to_datetime("2012-06-15")]

	## Aggregate the data
	## 2006
	names = ["punjab","sindh","khyber pakhtoon","balochistan"]
	dhs2006 = pd.Series([0.651, 0.507, 0.566, 0.463],
						index=["asia:pakistan:"+name for name in names])

	## 2012
	names = ["punjab","sindh","khyber pakhtoon","balochistan","islamabad","gilgit baltistan"]
	dhs2012 = pd.Series([0.70, 0.446, 0.578, 0.373, 0.852, 0.51],
						index=["asia:pakistan:"+name for name in names])

	return [dhs2006, dhs2012], dhs_years

if __name__ == "__main__":

	## Get the shapefile and create a sf object
	shp = "..\\_data\\Shapefile\\Pakistan.shp"
	sf = shapefile.Reader(shp) 
	kashmir = shapefile.Reader("..\\_data\\Shapefile\\polis_Kashmir.shp")
	kashmir = PathFromShape(kashmir.shapes()[0])

	## Get the DHS summary data
	series, years = DHSSummary()

	## Interpolate appropriately (just in time here, no space.)
	df = SpaceTimeDf(series,years)
	time = pd.DatetimeIndex(start="12/31/2008",end="12/31/2020",freq="SM")
	df = TimeInterpolate(df,freq="SM",time=time)
	
	## Pickle the result
	df.to_pickle("..\\pickle_jar\\extrapolated_ri.pkl")

	## Choose one a time to plot
	plot_time = "2017-11-30"
	coverage = df[:,plot_time]

	## Finally plot it up.
	fig, axes = plt.subplots(figsize=(10,8))
	#axes.add_patch(patches.PathPatch(kashmir,facecolor="None",alpha=1,lw=1))
	PlotDfonSf(fig,axes,coverage,sf,admin_level={1},alpha=1.,
			   colorbar=True,cmap=matplotlib.cm.Oranges_r,missing_color="grey")
	PlotBorders(fig,axes,sf,admin_level={0,1})
	axes.set(aspect="equal",title="MCV coverage, 2012 DHS, Pakistan")
	axes.axis("off")
	plt.tight_layout()
	plt.savefig("..\\_plots\\ri_map.png")

	## And plot a timeseries to
	## see the interpolation being used.
	admin = "asia:pakistan:punjab"
	timeseries = df.loc[admin]
	fig, axes = plt.subplots(figsize=(12,6))
	axes.plot(timeseries,color="0.4")
	axes.set(title="RI coverage, "+admin[len("asia:pakistan:"):].title())
	plt.show()

	