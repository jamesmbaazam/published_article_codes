# MCMC code to estimate the probability of a coin giving head

graphics.off(); # Close all figures left open from previous sessions
rm(list=ls()); # Clear all variables to start fresh!
tic = proc.time() # start timer
setwd("/Users/samik/Google Drive/Current Matlab/Lorenzo") # setting working directory on Mac (edit as required)

tic <- proc.time() # start timer


# Data given
n_tosses = 10
n_heads = 7

# If you want to generate other data, uncomment the following:
# p_true = 0.8;
# n_heads = binornd(n_tosses,p_true);

# Prior (uniform) - don't need to specify anything...

# MCMC parameters
niters = 10000 # Number of iterations
thinning = 1 # Save only 1 iteration every "thinning"
nsamples = floor(niters/thinning); # Effective number of samples that will be saved
# "floor" truncates to integer immediately below, i.e. floor(3.7623) = 3
sd_proposal = 0.5; # This is the width of the normal distribution used for the proposal
# Smaller "sd_proposal" mean more steps accepted but slower movement through parameter space
# Larger "sd_proposal" mean lots of rejection but broader exploration (hopping around a lot)

# Initialisation:
pchain = rep(0,nsamples); # Vector column of "nsample" zeros
LLchain = rep(0,nsamples); # Ditto

pc = runif(1); # Current (c) value of parameter p, drawn from the prior 
##### Can I fix this or should I draw it from prior?
LLc = log(dbinom(n_heads,n_tosses,pc)); # Current value of the log-likelihood (LL), i.e. calculated in current parameter value
# because the binomial coefficient is constant (i.e. indpendent of pc), we could also directly calculate this analytically as:
# LLc = n_heads * pc + (n_tosses-n_heads)*(1-pc)

accept_count = 0; # This is to count how many times I accept

for (s in 1:nsamples){ # start loop
  pchain[s] = pc; # Store current values of parameter
  LLchain[s] = LLc; # Store current value of the log-likelihood
  for (t in 1:thinning){
    pp = rnorm(1, pc, sd_proposal); # Proposed (p) new value for the parameter
    LLp = log(dbinom(n_heads,n_tosses,pp)); # Compute the new log-likelihood
    logalpha = LLp - LLc; # Compute the log of the probability of acceptance
    # logalpha has a particularly simple form because: 
    # 1) The prior is uniform (values cancel out), and 
    # 2) the proposal distribution is symmetric (values cancel out)
    if (logalpha >= 0 & !is.nan(logalpha)) { # i.e. if alpha >= 1 and is not NaN, accept for sure
      pc = pp # Save proposed value as new value
      LLc = LLp # Save log-likelihood of newly accepted value
      accept_count = accept_count + 1 # Record the acceptance event
    }
    else # if alpha < 1
      r = runif(1) # Generate a random number to decide whether to accept the new value or not
    if (log(r) < logalpha  & !is.nan(logalpha)) { # If r < alpha (and alpha is not NaN), accept (if in doubt, < is correct, because if alpha = 1 I would always accept)
      pc = pp # Save proposed value as new value
      LLc = LLp # Save log-likelihood of newly accepted value
      accept_count = accept_count + 1 # % Record the acceptance event
    }
  }
  cat(paste('Loop ', s, ' of ', nsamples, '\n'))
}

# PLOTTING #########

# To run only this section in R, click  anywhere below the "####" just
# above here and press "Code -> Run Region -> Run Code Section" (or, on Mac, alt+cmd+T)

# Parameters for post-MCMC analysis of the data
burnin = 1; # Define length of burnin period, i.e. how many samples to discard
nbins = 100; # Numer of bins for the histogram

# Trace plots
par(mfrow=c(2,1)) # Divide figure in 2 parts (2 rows, 1 column)
plot(burnin:nsamples, LLchain[burnin:nsamples],typ = "l", col = "black", lwd = 2,
     xlab = "Samples", ylab = "LL") # Plot log-likelihood (excluding burnin) against sample number
# xlab - write label to x-axis
# ylab - write label to y-axis
plot(burnin:nsamples, pchain[burnin:nsamples], typ = "l", col = "blue", lwd = 2,
     xlab = "Samples", ylab = "p") # Plot sampled values of p (excluding burnin)

# Autocorrelation
par(mfrow=c(1,1)) # Create a new figure, only one window
acf(pchain[burnin:nsamples],min(nsamples-1,100)); # Draw the autocorrelation plot
# This requires the function acf (built-in to R)

# Posterior distribution
par(mfrow=c(1,1)) # Create a new figure, only one window
hist(pchain[burnin:nsamples],nbins, freq=FALSE, main=NULL) # Plot the histogram of p 
# (excluding burnin), with "nbins" bins, and normalise it (rather than counting total numbers in each bin)
# If you want to zoom the x-axis only where the histogram is, uncomment the
# following 2 lines and comment the third one
# xl = par("usr"); xl = xl[1:2] # Get the limits for x-axis set automatically by matlab 
# x = seq(xl[1],xl[2],0.01); # Create a vector of value for the x axis, with step 0.01
x = seq(0, 1, 0.01);
abline(h = 1, lwd = 2, col = "red"); # Plot the prior, i.e. horizontal line at 1 - abline adds horizontal (h) or vertical (v) lines
# use "points" or "lines" instead of "plot" to add additional plots
lines(x, dbeta(x,1+n_heads,1+n_tosses-n_heads), typ="l", col="black"); # Plot the exact value of the posterior (which we know only because this is a simple case!)
ux <- unique(pchain[burnin:nsamples])
md = ux[which.max(tabulate(match(pchain[burnin:nsamples], ux)))] # calculate mode
title(paste('Mean = ', round(1e5*mean(pchain[burnin:nsamples]))/1e5, 
            '; mode = ', round(1e5*md)/1e5,
            '; sd = ', round(1e5*sd(pchain[burnin:nsamples]))/1e5,
            '; acceptance rate = ', round(1e4*100*accept_count/niters)/1e4, '%'))
# Write title, "mean" calculates mean, R doesn't have built-in "mode" so I did it manually,
# and "sd" to calculate standard deviation of the sample of p,
# and calculating the acceptance probability (written to 2 decimal places)

toc = proc.time() - tic # time taken to run it

toc

