""" tsir.py

A time series susepctible-infected-recovered model in python, which accounts for long term (linear) trends in 
susceptibles and for supplemental immunization campaigns (SIAs). This is based on the Grenfell 2000 work, and uses
generalized linear regression to calibrate the model."""

## For the linear algebra and
## data processing.
import numpy as np
import pandas as pd

## For interaction with the user
import warnings

################################################################################################
## TSIR model class
################################################################################################
class TSIR(object):

	""" Main TSIR model object, which performs maximum likelihood estimation given data, parametric 
	bootstrap resampling, and autogression for seasonality. """ 

	def __init__(self, dataframe, mcv1_take=0.9, mcv2_take=0.99):

		""" The object takes data as a dataframe, unpacks and processes it. The dataframe is expected to
		have the form

		df.columns = ['cases','population','birth_rate'] + optional (i.e. mcv1, mcv2, sia) 

		NB: All fractions, i.e. mcv and sia, are expected in decimal form, and there's currently no
		handling of nan's. Also birthrate must be in per time per 1000. people"""

		## Make sure the requirements are here
		assert ({"population","cases","birth_rate"} < set(dataframe.columns)), \
			   "Population, cases, and birth_rate are required - at least one is missing."

		## Keep the dataframe and the mcv takes
		## for reference purposes.
		self.df = dataframe
		self.mcv1_take = mcv1_take
		self.mcv2_take = mcv2_take

		## Certain methods are only available once 
		## the MLE has been computed. This is just to
		## keep track.
		self._mle_flag = False

		## Unpack the dataframe
		self.pop = dataframe["population"].values
		self.cases = dataframe["cases"].values
		self.birth_rate = dataframe["birth_rate"].values
		self.time = np.arange(len(self.pop))
		self.n_steps = len(self.time)

		## Create mcv1 and mcv2 columns from the dataframe,
		## otherwise set them to zero.
		try:
			assert (dataframe["mcv1"].max() <= 1.), "mcv1 needs to be from 0 to 1"
			self.mcv1 = dataframe["mcv1"].values
		except KeyError:
			self.mcv1 = np.zeros(self.pop.shape)
		try:
			assert (dataframe["mcv2"].max() <= 1.), "mcv2 needs to be from 0 to 1"
			self.mcv2 = dataframe["mcv2"].values
		except KeyError:
			self.mcv2 = np.zeros(self.pop.shape)

		## Create SIA column, which gets set to none
		## if not provided.
		try:
			sia_max = dataframe["sia"].max()
			assert (sia_max <= 1.), "sia coverage needs to be from 0 to 1"
			if sia_max == 0.:
				self.sia = None
				self.sia_times = None
			else:
				self.sia = dataframe["sia"].values
				self.sia_times = np.nonzero(self.sia != 0)[0]
		except KeyError:
			self.sia = None
			self.sia_times = None

		## Finally, compute the RI adjusted births
		births = self.birth_rate*(self.pop/1000.)
		self.adj_births = births*(1.-self.mcv1_take*self.mcv1*(1.-self.mcv2)-self.mcv2_take*self.mcv1*self.mcv2)
		
	def mle(self,weighted=True,verbose=False,detrended=True):

		""" Compute the MLE reporting rate, susceptible reconstruction, and infected reconstruction. This 
		is the outside facing method that runs 'private' methods depending on how the model has been 
		initialized. """

		## Trend fitting experiments failed, so now we assume
		## detrended susceptibles.
		if not detrended:
			raise NotImplementedError("Detrending is assumed now.")

		## We also force weigthed regression
		if not weighted:
			raise NotImplementedError("Susceptible reconstruction must be weighted.")

		if self.sia is None:
			self._mle_noSIA(verbose)
		else:
			self._mle_SIA(verbose)

		## internal boolean to tell some methods that 
		## the mle has been computed
		self._mle_flag = True

	def _mle_noSIA(self,verbose):

		""" No SIA version of MLE. Even though the SIA version would work, not computing the SIA
		matrix makes this faster. """

		## Construct the output of the model, in this case
		## the births+1
		self.output = np.cumsum(self.adj_births+1.)

		## And then create the feature matrix. Without SIA's, this 
		## is simply the cases plus 1.
		self.features = np.cumsum(self.cases+1.).reshape(-1,1)

		## Construct the weights. Weights are based on the variance
		## of [I_t | C_t, p].
		self.weights = 1./np.sqrt(self.cases + 1.)

		## Compute the MLE
		beta, beta_var, residual = WeightedLeastSquares(self.features,self.output,self.weights,verbose)

		## Store fit results
		self.beta = beta
		self.beta_var = beta_var
		self.Z = residual
		self.reporting_rate = 1./beta[0]

	def _mle_SIA(self,verbose):

		""" SIA version of the MLE, which constructs the SIA adjustment matrices and computes the
		fits accordingly. """

		## Construct the SIA matrix, a (t x t) matrix
		## with 1-sia coverage cumulative products as entries.
		self.M_sia = np.zeros((self.n_steps-1,self.n_steps-1))
		sia_c = 1. - self.sia
		for i in range(self.n_steps-1):
			self.M_sia[i:,i] = np.cumprod(sia_c[i:-1])

		## Now construct the 3 adjustment matrices. 	
		## D (for the intercept) and A (for the births and cases)
		## For D, we need a column of zeros and then the sia block.
		D = np.zeros((self.n_steps,self.n_steps))
		D[1:,1:] = self.M_sia

		## A has a repeated M_sia column as the first column, 
		## with a 1 in the extra entry. (See the notes for details).
		A = np.eye(self.n_steps)
		A[1:,0] = self.M_sia[:,0]
		A[1:,1:] = self.M_sia

		## Now we create the output vector
		self.output = np.dot(A,self.adj_births + 1.)

		## And the feature matrix. x0 is the sia-adjusted cumulative cases,
		## which corresponds to the reporting rate and x1
		## corresponds to the intercept (which is identifiable only due to
		## sia's).
		x0 = np.dot(A,(self.cases + 1.).reshape(-1,1))
		x1 = np.zeros((self.n_steps,1)) 
		x1[1:,0] = self.sia[:-1]/sia_c[:-1]
		x1 = np.dot(D,x1)

		## Since we're using the detrended method, S_t = S_bar + Z_t where
		## S_bar is a constant.
		self.features = np.concatenate([x0,x1],axis=1)

		## Construct the weights. weights are based on the Bayesian approach's
		## variance in the I_t | C_t distribution
		self.weights = 1./np.sqrt(self.cases + 1.)

		## Compute the MLE
		beta, beta_var, residual = WeightedLeastSquares(self.features,self.output,self.weights,verbose)

		## Store fit results
		self.beta = beta
		self.beta_var = beta_var
		self.Z = residual
		self.reporting_rate = 1./beta[0]

	def transmission_regression(self,periodicity=24):

		""" Regression of the transmission equation under the assumption that 
			
						I_t = beta_{t%periodicity}*(S_{t-1})*(I_{t-1}^alpha)*eps_t

			where epsilon is log-normally distributed and S_t = S_bar + Z_t. In principle, for
			S_bar = f(t), we could use the MLE of S_bar, but since that's not stable yet, it's assumed
			to be a constant fit parameter here (i.e. detrended MLE is assumed true). """

		assert self._mle_flag, "Need to have a detrended MLE to compute the transmission regression!"

		## Infered infecteds
		I_inferred = self.beta[0]*(self.cases+1.) - 1.

		## Create feature matrix for the auto regression,
		## as well as the step ahead output
		self.t_features = np.zeros((self.n_steps-1,periodicity+2))
		self.t_output = np.zeros((self.n_steps-1,))
		for i in range(self.n_steps-1):

			## Time of year
			time_in_period = (i+1) % periodicity
			self.t_features[i,time_in_period] = 1.

			## Infecteds in this time and the following
			## time.
			self.t_features[i,periodicity] = np.log(I_inferred[i])
			self.t_output[i] = np.log(I_inferred[i+1])

			## Susceptible residual
			self.t_features[i,periodicity+1] = self.Z[i]

		## Compute the regression
		params, params_var, residual = WeightedLeastSquares(self.t_features,self.t_output,standardize=False)

		## Use the residual to estimate Std dev of the data
		## in log space.
		RSS = (residual)**2
		var = RSS.sum(axis=0)/(self.t_features.shape[0] - self.t_features.shape[1])
		#mean_logE = np.sqrt(np.exp(var))
		std_logE = np.sqrt(var)#np.sqrt(np.exp(var)*(np.exp(var) - 1.))

		## Store the relevant pieces
		self.t_params = params
		self.t_var = params_var
		self.S_bar = 1./params[periodicity+1]
		self.t_beta = np.exp(params[:periodicity])/self.S_bar
		self.alpha = params[periodicity]
		self.periodicity = periodicity
		self.std_logE = std_logE


################################################################################################
## Linear regression stuff
def WeightedLeastSquares(X,Y,weights=None,verbose=False,standardize=True):

	""" Weighted LS, reduces to OLS when weights is None. This implementation computes
	the estimator and covariance matrix based on sample variance. For TSIR however, the

	NB: x is assumed to be an array with shape = (num_data points, num_features) """

	## Get the dimensions
	num_data_points = X.shape[0]
	try:
		num_features = X.shape[1]
	except:
		num_features = 1
		X = X.reshape((num_data_points,1))

	## Initialize weight matrix
	if weights is None:
		W = np.eye(num_data_points)
	else:
		W = np.diag(weights)

	## Standardize the inputs and outputs to help with
	## stability of the matrix inversion. This is needed because
	## cumulative cases and births both get very large.
	if standardize:
		muY = Y.mean()
		sigY = Y.std()
		muX = X.mean(axis=0)
		sigX = X.std(axis=0)
		X = (X-muX)/sigX
		Y = (Y-muY)/sigY

	## Compute the required matrix inversion
	## i.e. inv(x.T*w*x), which comes from minimizing
	## the residual sum of squares (RSS) and solving for
	## the optimum coefficients. See eq. 3.6 in EST
	xTwx_inv = np.linalg.inv(np.dot(X.T,np.dot(W,X)))

	## Now use that matrix to compute the optimum coefficients
	## and their uncertainty.
	beta_hat = np.dot(xTwx_inv,np.dot(X.T,np.dot(W,Y)))

	## Compute the estimated variance in the data points
	residual = Y - np.dot(X,beta_hat)
	RSS = (residual)**2
	var = RSS.sum(axis=0)/(num_data_points - num_features)

	## Then the uncertainty (covariance matrix) is simply a 
	## reapplication of the inv(x.T*x):
	beta_var = var*xTwx_inv

	## Reshape the outputs
	beta_hat = beta_hat.reshape((num_features,))

	## Rescale back to old values
	if standardize:
		X = sigX*X + muX
		Y = sigY*Y + muY
		beta_hat = beta_hat*(sigY/sigX)
		sig = np.diag(sigY/sigX)
		beta_var = np.dot(sig,np.dot(beta_var,sig))
		residual = sigY*residual + muY - np.dot(muX,beta_hat)

	## Print summary if needed
	if verbose:
		for i in range(num_features):
			output = (i,beta_hat[i],2.*np.sqrt(beta_var[i,i]))
			print("Feature %i: coeff = %.4f +/- %.3f." % output)

	return beta_hat, beta_var, residual