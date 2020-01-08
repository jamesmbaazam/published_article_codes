""" linelist_tools.py

Functions to manipulate line list data."""

## Standard stuff
import numpy as np
import pandas as pd

## For string matching
from fuzzywuzzy import process

####################################################################################################
## String matching
####################################################################################################
def InteractiveStandardizeNames(dot_names, df, match_table={}, drop_col=False, tol=99.,limit=50,
								limit_score=80.,preprocess = lambda s: s, dropna=False):

	"""Function to correct df.dot_name to as best as possible match the list of provided dot_names.
	This function uses fuzzy wuzzy's fuzzy string matching and assigns the best possible match. If 
	drop = True, the matched names overwrite the current dot names.

	match_table is a dictionary which maps provided names to dot names. It's an optional input if
	there's some match that you want to enforce. Otherwise, the match table is determined via fuzzy
	wuzzy.

	In cases where the match score is less than tolerance, limit number of options are computed and
	all with score > limit_score are offered to the user to resolve. 

	preprocess is a function which takes every df.dot_name and performs string manipulations on it. This is
	if you notice persistent issues in all strings. It's initialized to do nothing. For example, in the Nigeria 
	dataset, preprocess = lambda s: s.replace("_", " ") is a useful one.

	NB: This will no longer work if the input isn't from GetRawData (i.e. if it has a multi index for
	space and time). Need to fix this."""

	assert ~isinstance(df.index, pd.core.index.MultiIndex), \
			"Can't interact with a Multiindex data frame, needs to be the raw data."

	## Storage for the matched names columns
	matched_names = []

	## Loop through provided names 
	## and compute matches as needed.
	for provided_name in df.dot_name:

		## Check to see if we have this match already
		## or if it's been provided.
		match = match_table.get(provided_name, None)

		## If we don't have it, use fuzzy matching to find it.
		if match is None:
			
			## Get the best match
			match, score = process.extractOne(preprocess(provided_name), dot_names)
			
			## For scores below tolerance, we print limit options and 
			## ask the user to resolve the discrepency.
			if score <= tol:	
				print("\n"+provided_name+" has a match that scores under tolerance.")
				print("Here are potential options:")
				
				## Compute and print the options.
				matches = process.extract(preprocess(provided_name), dot_names, limit=limit)
				for i, m in enumerate(matches):
					if m[1] < limit_score:
						break
					print("{}: {}".format(i,m))

				## Get user match. If it can be cast as an interger,
				## it's interpreted as an index in matches. If not, it's
				## taken as is. 
				match = input("Which name do you want matched with "+provided_name+" (if none, type nan)? ")
				try: 
					match = matches[int(match)][0]
				except ValueError:
					match = match

				## Verify what was selected to prevent mistakes.
				print("You selected "+match)

			## Save the match.
			match_table[provided_name] = match
		
		## Store the match for the df
		matched_names.append(match)

	## Check performance by testing the
	## number of unused names.
	used = set(match_table.values())
	unused = set(dot_names) - used

	## Overwrite dot_names if required.
	if drop_col:
		df["dot_name"] = matched_names
	else:
		df["matched_name"] = matched_names

	## If we want to drop nans
	if dropna:
		df = df.replace("nan",np.nan).dropna()
	
	return df, match_table, unused

####################################################################################################
## Other useful tools
####################################################################################################
def Accumulate(df,freq,admin_col="dot_name",name="cases",t_index=None):

	## Create uniform time index for each admin unit
	if t_index is None:
		t_index = pd.DatetimeIndex(start=df["time"].min(),end=df["time"].max(),freq=freq)
	labels = t_index[:-1]

	## Group by dot name 
	data = {}
	for name, subframe in df.groupby(admin_col):
		bins = pd.cut(subframe["time"],t_index,labels=labels)
		data[name] = bins.value_counts().reindex(labels).fillna(0.)

	## Create a data frame with the appropriate names	
	new_df = pd.concat(data.values(), keys=data.keys())
	new_df.index.rename(["admin","time"],inplace=True)
	return new_df


if __name__ == "__main__":

	print("See country specific directories for examples of usage, specifically the LineList.py files.")