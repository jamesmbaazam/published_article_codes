""" NightTimeLights.py

Scripts and processing to create monthly night time lights timeseries using rasters from 
https://ngdc.noaa.gov/eog/viirs/download_dnb_composites.html

This is to quanitify population fluctuations in major cities (Karachi, Lahore). """
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

## Risk map tools
from riskmap3.map_maker import *
plt.rcParams["font.size"] = 24.

## For getting file names
from os import listdir

if __name__ == "__main__":

	## Get the shapefile and create a sf object
	shp = "..\\_data\\Shapefile\\Pakistan_majorcities.shp"
	sf = shapefile.Reader(shp)

	## And the country shapeifle
	admin0 = shapefile.Reader("..\\_data\\Shapefile\\polis_PAK_country_revised.shp")
	bbox = admin0.shapes()[0].bbox
	points = np.array(admin0.shapes()[0].points)

	## Get the raster
	root = "D:\\nthakkar\\PAK_lights\\"
	fnames = listdir(root)

	## Loop over rasters and then shapes and compute total brightness
	brightness = []
	plot = False

	## Get the date from the fname
	for fname in fnames:

		## Extract the date
		date = pd.to_datetime(fname[fname.find("npp_")+4:fname.find("-")],format="%Y%m%d")
	
		for sr in sf.shapeRecords():

			## Get the dotname
			dot_name = sr.record[1].lower()
			
			## Get the associated part of the raster
			bbox = sr.shape.bbox
			extent = (bbox[0],bbox[2],bbox[1],bbox[3])
			raster = georaster.SingleBandRaster(root+fname,load_data=extent,latlon=True)

			## Compute total brightness and size of the
			## bounding box
			num_pixels = int(np.prod(raster.r.shape))
			total_brightness = raster.r.sum()

			## Append the result
			entry = (dot_name,date,total_brightness,num_pixels)
			brightness.append(entry)

			if plot and dot_name == "asia:pakistan:punjab:lahore":
				plt.imshow(raster.r,extent=raster.extent)
				points = np.array(sr.shape.points)
				plt.plot(points[:,0],points[:,1],color="k",lw=2)
				plt.show()

	## Create a dataframe
	df = pd.DataFrame(brightness,columns=["dot_name","time","total","num_pixels"])
	df = df.set_index(["dot_name","time"]).sort_index(level=[0,1])
	print(df)

	## Pickle the result
	df.to_pickle("..\\pickle_jar\\nighttimelights.pkl")
	
	## Groupby month of the year for brightness traces
	df = df.groupby(["dot_name",lambda s: s[1].month]).mean()
	df["per_pixel"] = df["total"]/df["num_pixels"]

	## Find the total mean
	mean = df.groupby(level=1)["per_pixel"].mean()
	std = df.groupby(level=1)["per_pixel"].std()

	## Make a plot
	fig, axes = plt.subplots(figsize=(12,8))
	for dn, subframe in df.groupby("dot_name"):
		axes.plot(subframe.loc[dn,"per_pixel"],color="grey",alpha=0.5,lw=1)
	#axes.fill_between(mean.index,mean-std,mean+std,color="C0",alpha=0.2)
	axes.plot(mean,color="k",lw=3)
	axes.set(xlabel="Month of the year",ylabel="Brightness per pixel")
	plt.tight_layout()
	plt.show()


