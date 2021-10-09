""" vis_tools.py

Functions to make plots that summarize and help interpret the data."""

## Standard stuff
import numpy as np
import pandas as pd

## For plotting
import matplotlib.pyplot as plt

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
plt.rcParams["mathtext.fontset"] = "dejavuserif"

def DataSummary(df,admin_name=None):

	""" Create a plot of all the time series at a specific admin name 
	within the df. """

	## What's the number of columns, that'll
	## set the number of rows in the plot.
	columns = df.columns
	num_rows = len(columns)

	## Create the plot
	figsize = (12,5+num_rows)
	fig, axes = plt.subplots(num_rows,1,sharex=True,figsize=figsize)
	for i,column in enumerate(columns):
		color = "C"+str(i % 9)
		axes[i].plot(df[column],label=column,color=color)
		axes[i].legend()
		axes[i].set(ylabel=column)
	axes[-1].set(xlabel="Time")
	if admin_name:
		axes[0].set(title="Summary for "+admin_name)

	return fig, axes