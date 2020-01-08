""" timeseries_tools.py

Functions to facilitate time series analysis by taking lists of pandas series with
spactial indexes and combining/interpolating/concatenating into multi-index space time 
data frames."""

## These functions mostly just wrap pandas timeseries
## functionality. So we don't need much.
import pandas as pd

def SpaceTimeDf(series,years,level_names=["admin","time"]):

	"""Function to make a multiindex df from list of series at specific years (also list)."""

	assert all(isinstance(y,pd.Timestamp) for y in years), "Years needs to be a pandas object. Use pd.to_datetime()"

	## Create a data frame with the appropriate names
	df = pd.concat(series, keys=years).swaplevel(i=0,j=1).sort_index(level=0)
	df.index.rename(level_names,inplace=True)
	return df

def EmptySpaceTimeDf(names,time_index,dtype=float,level_names=["admin","time"],val=0.):

	""" Function to create an initialized Space time df. """

	empty = {name:pd.Series(len(time_index)*[val],index=time_index,dtype=dtype) for name in names}
	empty = pd.concat(empty.values(),keys=empty.keys())
	empty.index.rename(level_names,inplace=True)
	return empty

def TimeInterpolate(df,freq,time=None,groupby="admin",level_names=["admin","time"]):

	"""Function to, by admin level, resample the df along given times or at a specfic frequency. As is
	this constructs a totally new df with the same structure. This isn't ideal, there should be a way to
	work in place within the given df."""

	new_df = {}
	for admin, sf in df.groupby(groupby):
		timeseries = sf[admin,:].resample(freq).mean().interpolate()
		if time is not None:
			timeseries = timeseries.reindex(time,method="nearest").fillna(method="ffill").fillna(method="bfill")
		new_df[admin] = timeseries

	new_df = pd.concat(new_df.values(),keys=new_df.keys())
	new_df.index.rename(level_names,inplace=True)
	return new_df