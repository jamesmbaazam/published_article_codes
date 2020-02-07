# campaign_timing
A time series model for forecasting measles vaccination campaign impact.

This is a repository with Python code pertaining to the paper *Decreasing measles burden by optimizing campaign timing* by Niket Thakkar, Syed Saqlain Ahmad Gilani, Quamrul Hasan and Kevin McCarthy, 2019.

The folder `riskmap3` contains the basic model class and a set of data analysis functions. Those objects are used for data processing and modelling in the scripts within the `Pakistan` directory.

Although all relevant code is provided, data from Pakistan is highly sensitive cannot be shared on Git. As a result, all code makes reference to data within the `_data` directory and serialized pandas dataframes within the `pickle_jar` that are not provided. For more information, please see the open source manuscript.
