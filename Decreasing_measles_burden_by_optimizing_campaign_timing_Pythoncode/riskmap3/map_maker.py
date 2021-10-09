""" map_maker.py

Functions to facilitate plotting pandas dataframes onto shapefiles."""

## Standards
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## For shape files
import shapefile

## For plotting shapefile polygons
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

## Custom global matplotlib parameters
## see http://matplotlib.org/users/customizing.html for details.
plt.rcParams["font.size"] = 18.
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["font.serif"] = "Garamond"
plt.rcParams["xtick.labelsize"] = "medium"
plt.rcParams["ytick.labelsize"] = "medium"
plt.rcParams["legend.fontsize"] = "medium"
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["axes.formatter.use_mathtext"] = True
plt.rcParams["mathtext.fontset"] = "cm"

##############################################################################################################
## Plotting functions
##############################################################################################################
def PlotDfonSf(fig,axes,df,sf,
			admin_level={0,1,2,3},cmap=matplotlib.cm.viridis,colorbar=True,linecolors={0:'0',1:'0.12',2:'0.25',3:'0.5'},
			alpha=1.,vlim=None,missing_color="white"):

	'''Plots a pandas df on a shapefile using ploygons in matplotlib.

	df: simple pandas series where index corresponds to records in the shapefile and values are to be plotted
	sf: full shapefile object from pyshp'''

	## Set the color bar limits based
	## on the df if none is provided.
	if vlim is None:
		vlim = (df.min(),df.max())

	## Storage to create a patch collection
	patches = []
	colors = []

	for ShapeRecord in sf.iterShapeRecords():

		## Unpack the Shape+Record object
		shape = ShapeRecord.shape
		record = ShapeRecord.record

		## Figure out the level of the shape based on
		## the number of colons in the name string
		name_str = record[1].lower()
		this_level = name_str.count(":")-1

		## Pass over the shape if it's not in the 
		## desired admin levels.
		if this_level not in admin_level:
			continue

		## Figure out if the shape is broken into parts
		## and add a -1 to be able to index appropriates
		parts = shape.parts
		parts.append(-1) 

		## Loop through each part and create seperate polygons
		## for each.
		points = shape.points
		points.append(points[0])
		points = np.array(points)
		for i,j in zip(parts[:-1],parts[1:]):

			## Create a polygon based shape.points
			polygon = Polygon(points[i:j,:])
			patches.append(polygon)

			## Get the face color from the dataframe. If not available,
			## use -np.inf and later specify a set_under value for the cmap.
			try:
				colors.append(df[name_str])
			except:
				colors.append(-np.inf)

	## Create matplotlib patch collection
	## with specified colorbar
	cmap.set_under(missing_color)
	collection = PatchCollection(patches, cmap=cmap, alpha=alpha)
	collection.set_array(np.array(colors))
	collection.set_clim(vlim)

	## Finally plot
	axes.add_collection(collection)
	if colorbar:
		fig.colorbar(collection, ax=axes, fraction=0.035)

	## And set the limits via matplotlib autoscale
	axes.autoscale_view()

def PlotBorders(fig,axes,sf,admin_level={0,1,2,3},colors={0:'0',1:'0.12',2:'0.25',3:'0.5'}):

	"""Function to make a simple plot of the boundaries in the shape file.

	sf: the shapefile object to plot
	admin_level: set of admin levels to plot (0: country, 1: state, 2: LGA, etc.)
	colors: dictionary that maps admin-level to color. """

	for ShapeRecord in sf.iterShapeRecords():

		## Unpack the Shape+Record object
		shape = ShapeRecord.shape
		record = ShapeRecord.record

		## Figure out the level of the shape based on
		## the number of colons in the name string
		name_str = record[1]
		this_level = name_str.count(":")-1

		## Pass over the shape if it's not in the 
		## desired admin levels.
		if this_level not in admin_level:
			continue
	
		## Figure out if the shape is broken into parts
		## and add a -1 to be able to index appropriates
		parts = shape.parts
		parts.append(-1) 

		## Create a np array of the points list to facilitate matplotlib.
		## We append the first point to make the boundary a full loop and to be
		## able to use -1 as the end below.
		points = shape.points
		points.append(points[0])
		points = np.array(points)

		## Plot by individual parts
		for i,j in zip(parts[:-1],parts[1:]):
			axes.plot(points[i:j,0],points[i:j,1],c=colors.get(this_level,"0"),lw=max(1.5 - 0.5*this_level,0.5))


if __name__ == "__main__":

	### Example usage of the above function.

	## Get the shapefile and create a sf object
	shp = "..\\_data\\Shapefile\\Nigeria.shp"
	sf = shapefile.Reader(shp)

	## Create a dataframe to be plotted on the map
	names = [record[1] for record in sf.records() if record[1].count(":") == 3]
	df = pd.Series({n.lower():i for n,i in zip(names,range(len(names)))})

	## Make the plot
	fig, axes = plt.subplots()
	PlotDfonSf(fig,axes,df,sf,colorbar=True,admin_level={2},alpha=0.45)
	PlotBorders(fig,axes,sf,admin_level={0,1,2})
	plt.show()
