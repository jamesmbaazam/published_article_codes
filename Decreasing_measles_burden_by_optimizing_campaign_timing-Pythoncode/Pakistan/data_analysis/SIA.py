""" SIA.py

Coverage and dates for the SIA call. This is done more or less "by hand" since SIA calendars
are in wild formats."""
import sys
sys.path.insert(0,"..\\..\\")

import numpy as np
import pandas as pd
import shapefile
from riskmap3.data_process.shapefile_tools import NamesFromSf

def SIASummary():

	""" Data from Measles SIA 2007_2015_administrative report.xls, 
	processed by hand, not automatically. 

	NB: efficacy is the amount of the base coverage that SIA counts for (based on how much
	pop was targeted). Numbers here came from Kurt."""

	## For date time and dot name conversion.
	to_dt = lambda t: pd.to_datetime(t,format="%d-%m-%Y")
	to_dn = lambda ps: ["asia:pakistan:"+p for p in ps]
	base_coverage = 0.4
	future_coverage = 1.

	## Storage for the SIA series
	sias = []
	
	## 2010 SIA
	provinces = ["punjab","sindh","khyber pakhtoon","balochistan"]
	start_dates = pd.Series(["20-09-2010","25-10-2010","20-09-2010","20-09-2010"],
							index=to_dn(provinces),name="start_date").apply(to_dt)
	end_dates = pd.Series(["02-10-2010","05-11-2010","02-10-2010","02-10-2010"],
							index=to_dn(provinces),name="end_date").apply(to_dt)
	efficacy = pd.Series([0.18,0.39,0.6,0.1],index=to_dn(provinces),name="coverage")
	sia2010 = pd.concat([start_dates,end_dates,base_coverage*efficacy],axis=1)
	sias.append(sia2010)

	## 2011 SIA
	provinces = ["punjab","sindh","khyber pakhtoon","balochistan","sindh","gilgit baltistan","ajk"]
	start_dates = pd.Series(["17-01-2011","05-01-2011","10-01-2011","05-01-2011","18-02-2011","04-07-2011","21-11-2011"],
							index=to_dn(provinces),name="start_date").apply(to_dt)
	end_dates = pd.Series(["29-01-2011","17-01-2011","22-01-2011","17-01-2011","01-03-2011","09-07-2011","26-11-2011"],
							index=to_dn(provinces),name="end_date").apply(to_dt)
	efficacy = pd.Series([0.36,0.18,0.25,0.11,0.42,1.,1.],index=to_dn(provinces),name="coverage")
	sia2011 = pd.concat([start_dates,end_dates,base_coverage*efficacy],axis=1)
	sias.append(sia2011)

	## 2012 SIA (SHOULD THERE BE MORE HERE?) The dates are unclear, I'm using
	## the MSL SIA 2012 dates for balochistan for everywhere.
	provinces = ["balochistan"]#,"punjab","khyber pakhtoon","fata","islamabad"]
	start_dates = pd.Series(["24-12-2012"]*len(provinces),index=to_dn(provinces),name="start_date").apply(to_dt)
	end_dates = pd.Series(["04-01-2013"]*len(provinces),index=to_dn(provinces),name="end_date").apply(to_dt)
	efficacy = pd.Series([0.81],index=to_dn(provinces),name="coverage")
	sia2012 = pd.concat([start_dates,end_dates,base_coverage*efficacy],axis=1)
	sias.append(sia2012)

	## 2013 SIA
	## NOT SURE WHAT TO DO HERE...(these dates came from Kurt)
	provinces = ["punjab","sindh"]
	start_dates = pd.Series(["15-05-2013"]*len(provinces),index=to_dn(provinces),name="start_date").apply(to_dt)
	end_dates = pd.Series(["15-05-2013"]*len(provinces),index=to_dn(provinces),name="end_date").apply(to_dt)
	efficacy = pd.Series([0.66,0.66],index=to_dn(provinces),name="coverage")
	sia2013 = pd.concat([start_dates,end_dates,base_coverage*efficacy],axis=1)
	sias.append(sia2013)

	## 2014 SIA (what's the difference between CDA and islamabad? I'm gonna combine those SIAs into
	## one).
	provinces = ["punjab","sindh","khyber pakhtoon","balochistan","ajk","gilgit baltistan","fata","islamabad"]
	start_dates = pd.Series(["26-01-2015","19-05-2014","19-05-2014","13-04-2015","08-12-2014","18-05-2015",
							"20-08-2015","09-02-2015"],index=to_dn(provinces),name="start_date").apply(to_dt)
	end_dates = pd.Series(["09-02-2015","31-05-2014","31-05-2014","25-04-2015","20-12-2014","23-05-2015",
						   "31-08-2015","28-02-2015"],index=to_dn(provinces),name="end_date").apply(to_dt)
	efficacy = pd.Series([1.]*len(provinces),index=to_dn(provinces),name="coverage")
	sia2014 = pd.concat([start_dates,end_dates,base_coverage*efficacy],axis=1)
	sias.append(sia2014)

	## 2017 SIA
	provinces = ["sindh"]
	start_dates = pd.Series(["15-08-2017"],index=to_dn(provinces),name="start_date").apply(to_dt)
	end_dates = pd.Series(["15-08-2017"],index=to_dn(provinces),name="end_date").apply(to_dt)
	efficacy = pd.Series([0.15],index=to_dn(provinces),name="coverage")
	sia2017 = pd.concat([start_dates,end_dates,base_coverage*efficacy],axis=1)
	sias.append(sia2017)

	#################### Future SIAs??
	provinces = ["punjab","sindh","khyber pakhtoon","balochistan","ajk","gilgit baltistan","fata","islamabad"]

	## Each province is visted twice
	provinces = [p for pair in zip(provinces,provinces) for p in pair]

	## Start and end
	start_dates = pd.Series(["01-12-2019","17-12-2019"]*int(len(provinces)/2),index=to_dn(provinces),name="start_date").apply(to_dt)
	end_dates = pd.Series(["11-12-2019","27-12-2019"]*int(len(provinces)/2),index=to_dn(provinces),name="end_date").apply(to_dt)
	efficacy = pd.Series([0.5]*len(provinces),index=to_dn(provinces),name="coverage")
	future_sia = pd.concat([start_dates,end_dates,future_coverage*base_coverage*efficacy],axis=1)
	sias.append(future_sia)

	## Send back a full summary
	return pd.concat(sias,axis=0)

if __name__ == "__main__":

	## Get the shapefile
	sf = shapefile.Reader("..\\_data\\Shapefile\\Pakistan.shp")

	## Get the data
	sias = SIASummary()

	## Compute a midpoint for the sia time in the model.
	sias["mid_point"] = sias.start_date + 0.5*(sias.end_date - sias.start_date)
	sias["mid_point"] = sias.mid_point.dt.round("d")
	sias = sias.sort_values("mid_point")
	
	## Now create the appropriate spacetime-df
	data = {}
	time = pd.DatetimeIndex(start="12/31/2008",end="12/31/2020",freq="SM")
	for province, df in sias.groupby(level=0):
		labels = pd.cut(df.mid_point,time,labels=time[1:])
		binned = labels.value_counts().reindex(time[1:])
		binned[binned != 0.] = binned[binned != 0.]*df["coverage"].as_matrix()
		data[province] = binned
	
	df = pd.concat(data.values(),keys=data.keys())
	df.index.rename(["admin","time"],inplace=True)

	#df.to_pickle("..\\pickle_jar\\extrapolated_sia_12_19.pkl")
	print(df[df != 0].to_latex())
	






