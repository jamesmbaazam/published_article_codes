""" Using Hil's script to calculate birthrates from the dhs."""

## Orient the script's path
import sys
sys.path.insert(0,"..\\..\\")

import numpy as np
import pandas as pd

import shapefile
from riskmap3.map_maker import *
from riskmap3.data_process.shapefile_tools import *
from riskmap3.data_process.timeseries_tools import *

## Get the shapefile and create a sf object
shp = "..\\_data\\Shapefile\\Pakistan.shp"
sf = shapefile.Reader(shp) 

## Get the output from Hil's script
from_hil = pd.read_csv("..\\_data\\BirthRateDHS\\DAT_OUT2.csv")
birth_rates = from_hil["yearly_birthrate_per_1000"].as_matrix()[:-1]

## Convert to monthly
birth_rates = 1000.*((1. + birth_rates/1000.)**(1./24.) - 1.)

## Change the index to the sf names, etc.
names = ["punjab","sindh","khyber pakhtoon","balochistan","gilgit baltistan","islamabad"]
dhs2012 = pd.Series(birth_rates,index=["asia:pakistan:"+name for name in names])

## Format as a space-time df
df = SpaceTimeDf([dhs2012],[pd.to_datetime("2012-06-15")])
time = pd.DatetimeIndex(start="01/01/2009",end="12/31/2020",freq="SM")
df = TimeInterpolate(df,freq="SM",time=time)

## Pickle the results
df.to_pickle("..\\pickle_jar\\extrapolated_dhs_birth_rate.pkl")

## Choose a time to plot
plot_time = "2017-11-30"
coverage = df[:,plot_time]

## Finally plot it up.
fig, axes = plt.subplots(figsize=(10,8))
PlotDfonSf(fig,axes,coverage,sf,admin_level={1},alpha=0.9,
		   colorbar=True,missing_color="grey")
PlotBorders(fig,axes,sf,admin_level={0,1})
axes.set(title="Birth rate, "+plot_time)
axes.set(title="Population, "+plot_time,aspect="equal")
axes.axis("off")
plt.tight_layout()

## And plot the timeseries
admin = "asia:pakistan:punjab"
timeseries = df.loc[admin]
fig, axes = plt.subplots(figsize=(12,6))
axes.plot(timeseries,color="0.4")
axes.set(title="Birth rate, "+admin[len("asia:pakistan:"):].title())
plt.show()