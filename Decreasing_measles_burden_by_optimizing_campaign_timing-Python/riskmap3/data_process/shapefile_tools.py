""" shapefile_tools.py

Functions to facilitate mapping gps data onto shape files, assigning admin names, etc.

NOTE: Maybe I should make a shapefile object that has all this stuff as methods to prevent
repeated computing of paths, etc. This is a short-term development goal."""

## Standards
import numpy as np
import pandas as pd

## For shape files
import shapefile

## For determining whether things are within or
## outside shapes
from matplotlib.path import Path

##############################################################################################################
## Functions for shape file interactions
##############################################################################################################
def AdminLevelSubset(sf,admin_level):

	""" Function to extract a dictionary of shape objects corresponding to a particular admin
	level. Key is dot name and value is shape object. """

	num_colons = admin_level+1
	shapes = {sr.record[1].lower():sr.shape for sr in sf.iterShapeRecords() if sr.record[1].count(":") == num_colons}
	return shapes

def NamesFromSf(sf,admin_level):

	""" Function to retrieve all the names in shape file at a given admin level. """

	num_colons = admin_level + 1
	names = [r[1].lower() for r in sf.iterRecords() if r[1].count(":") == num_colons]
	return names

def PathFromShape(shape):

	"""Create a matplotlib path object from a shape object."""

	## Start by creating the appropriate codes,
	## telling the path object when the shape is disconnected.
	codes = np.array(len(shape.points)*[Path.LINETO])
	codes[shape.parts] = Path.MOVETO

	## Create the path using shape.points
	path = Path(shape.points,codes)

	return path

def PointsWithin(test_points,shape):

	""" Check whether points are within a shape. Shape is a shapefile object, test_points is an 
	array like object with shape = (num_points, 2) of (x,y) coordinates. Returns a bool array
	of the shape = (num_points,) 
	
	NB: There are some known issues with edge cases here. A point on the edge of a shape will not
	return True, and for now, there's no way to change this."""

	## Create path object
	path = PathFromShape(shape)

	## Find points within
	within = path.contains_points(test_points)

	return within

def ShapeCenter(shape):

	"""Return the center point of the shape. For now this is done by treating the shape as a collection of points
	and finding the point that minimizes the L2 distance between all points. That, thankfully, is just the mean of the
	points. I'm not sure this is the right thing to do though, since sometimes that point can be outside the shape, especially
	for L-shaped units, islands, etc."""

	points = np.array(shape.points)
	return np.mean(points,axis=0)

def AdminLevelFromGPS(df,sf,admin_level,verbose=False,
					  col_cluster="v001",col_x="LONGNUM",col_y="LATNUM",output_col=None):

	"""Function to create a column in dataframe df (which is assumed to have columns col_x, col_y, and col_cluster),
	that has the name of admin_level unit which contains the points. This function handles points outside any admin
	unit via minimum L2 distance."""

	## Get admin unit shapes, which we 
	## store as a list of tuples to be able to
	## manipulate the order it's looped through.
	shapes = AdminLevelSubset(sf,admin_level)
	list_of_shapes = list(shapes.items())

	## Create dictionary of paths
	paths = {}

	## And of clusters
	clusters = {}

	## distance function and dictionary of centers
	## to handle unclear cases.
	centers = {name:ShapeCenter(admin) for name,admin in shapes.items()}
	distance = lambda x1,x2: np.linalg.norm(x1-x2)
	via_distance = 0

	## Loop through clusters and assign admin unit
	for v001, sub_frame in df.groupby(col_cluster):

		## Extract GPS
		point = np.array([[sub_frame[col_x].iloc[0],sub_frame[col_y].iloc[0]]])

		## Loop through possible shapes
		for i, name_shape in enumerate(list_of_shapes):

			## Unpack the name, shape tuple
			name, admin = name_shape

			## Bounding box check, to avoid checking
			## ones that aren't close.
			bbox = admin.bbox
			within_bbox = (bbox[0] <= point[0,0] <= bbox[2]) & (bbox[1] <= point[0,1] <= bbox[3])
			if not within_bbox:
				continue

			## If it's within the bounding box,
			## we need a more careful check. We can us
			## matplotlib for that.
			path = paths.get(name,None)
			if path is None:
				path = PathFromShape(admin)
				paths[name] = path

			## Is it inside?
			within = path.contains_points(point)[0]

			## If so, we store it and we rearrange our
			## list of shapes so this shape is guessed first, since
			## often neighboring points are in the same LGA (but that 
			## in general depends on the structure of the input).
			if within:
				clusters[v001] = name
				if i != 0:
					list_of_shapes.insert(0,list_of_shapes.pop(i))
				break

		## If we get here, the cluster was jittered
		## outside of all units at admin_level (and into the water or
		## something). So we can use the centriods to 
		## determine the best assignment.
		if clusters.get(v001,None) is None:
			via_distance += 1
			distances = {name:distance(point,center) for name, center in centers.items()}
			closest = min(distances, key=distances.get)
			clusters[v001] = closest

	## Create admin unit column
	if output_col is None:
		output_col = "admin"+str(admin_level)
	df[output_col] = [clusters.get(v001,"NaN") for v001 in df[col_cluster]]

	## Print summary
	if verbose:
		total_num = len(shapes)
		num_used = len(df[output_col].value_counts())
		print("Total number of admin %i = %i." % (admin_level,total_num))
		print("Number of admin %i with clusters = %i." % (admin_level,num_used))
		print("%i clusters were assigned via distance." % via_distance)
		print("%i admin %i have no clusters in them." % (total_num - num_used, admin_level))

	return df

def GridShape(shape,n=50):

	"""Return a nxn grid in the shape's bounding box, subset such that
	points outside the box are neglected. This is a dumb implementation of this,
	and will fail for pathelogically long and skinny shapes. Still, it works well for
	states and countries."""

	## Create 1d grids
	bbox = shape.bbox
	x_grid = np.linspace(bbox[0],bbox[2],n)
	y_grid = np.linspace(bbox[1],bbox[3],n)

	## Mesh them
	X,Y = np.meshgrid(x_grid, y_grid)

	## Find points within the shape and outside the shape
	points = np.array([X.reshape(-1), Y.reshape(-1)]).T
	within = PointsWithin(points, shape)

	## Return points within
	return points[within]

def NearestNeighbors(sf,admin_level,fill_empty=False):

	"""Function to calculate a dictionary which maps admin units to names of neighbors."""
	
	## Get the admin level subset
	admins = AdminLevelSubset(sf,admin_level)

	## Calculate matplotlib paths
	paths = {name:PathFromShape(admin) for name, admin in admins.items()}

	## Assume for now that no islands exist. This will
	## be corrected as we compute neighbors.
	islands = []

	## Loop through paths and create a dictionary of neighbors
	neighbors = {}
	for n1, p1 in paths.items():
		
		## List of names of this admin unit's
		## neighbors.
		this_set = []
		
		## Loop through other paths. In principle, this
		## can be optimized since being a neighbor implies
		## having a neighbor, which we take advantage of. But
		## the speed gains are modest.
		for n2, p2 in paths.items():
			if n1 == n2:
				continue
			if n1 in neighbors.get(n2,[]):
				this_set.append(n2)
				continue
			if p1.intersects_path(p2):
				this_set.append(n2)

		## Store it
		if not this_set:
			islands.append(n1)
		neighbors[n1] = this_set

	## Correct islands if required, both by user
	## and by the shape file itself. In these cases, neighbors are
	## assumed to be the neighbors of the L2 distance nearest neighbor.
	if islands and fill_empty:

		## Start by computing admin-unit centers, and specifying
		## the distance function.
		centers = {name:ShapeCenter(admin) for name,admin in admins.items()}
		distance = lambda x1,x2: np.linalg.norm(x1-x2)

		## Now find min distances
		for island in islands:
			
			## Sort centers by distance
			this_center = centers[island]
			nearest = sorted(centers.keys(), key=lambda x: distance(this_center,centers.get(x)))

			## Find the first non-empty neighbor
			for name in nearest[1:]:
				if neighbors[name]:
					neighbors[island] = neighbors[name]+[name]
					print("Matching "+island+" with "+name+".")
					break

	return neighbors


if __name__ == "__main__":

	shp = "..\\..\\_data\\Shapefile\\Nigeria.shp"
	sf = shapefile.Reader(shp)
	#print(len(NamesFromSf(sf,2)))

	
	from datetime import datetime
	start = datetime.now()

	neighbors = NearestNeighbors(sf,admin_level=2,fill_empty=False)
	end = datetime.now()
	print(end-start)
	#print(neighbors)

	for n, lst in neighbors.items():
		if not lst:
			print(n)
	

