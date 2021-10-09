# MCMC code to estimate transmission and recovery rates from the final size
# of an SIR epidemic in a large population. The final size is measured by
# sampling "n_tested" individuals (with a perfect test) and detecting that 
# "n_pos" have been infected

# if the program complains about not having the function "mvrnorm", run these following two lines once
# install.packages("MASS") 
# library("MASS")

graphics.off(); # Close all figures left open from previous sessions
rm(list=ls()); # Clear all variables to start fresh!
tic = proc.time() # start timer
setwd("~/Google Drive/Current Matlab/Lorenzo") # correct folder
source("Gamma_meansd2shapescale.r") # function for gamma mean and sd
source("LL_EpidemicFinalSizeLargePop.r") # function for likelihood
source("scatterhist.r") # function for likelihood


# Data given
n_tested = 100;
n_pos = 80;

# Propose new parameters in block or one at a time?
block_proposal = TRUE; 
# if true, propose both parameters from a bivariate normal distribution
# if false, propose one parameter at a time

# Priors:
# For beta: flat, improper prior on [0,+infinity)
# For gamma: gamma-distributed with mean m_hy_gamma days and sd sd_hy_gamma days
m_hy_gamma = 4; # hyper-parameter for shape of Gamma distribution for gamma
sd_hy_gamma = 1; # hyper-parameter for scale of Gamma distribution of gamma
# R parameterises the gamma distribution in terms of shape and scale parameters
bla = Gamma_meansd2shapescale(m_hy_gamma,sd_hy_gamma); # Find corresponding shape and scale
bla = sapply(bla, function(x){as.numeric(x[1])})
shape = bla[1]
scale= bla[2]

# MCMC parameters
niters = 10000; # Number of iterations
thinning = 1; # Save only 1 iteration every "thinning"
nsamples = floor(niters/thinning); # Effective number of samples that will be saved
# "floor" truncates to integer immediately below, i.e. floor(3.7623) = 3
sd_proposals = c(2, 1); # This is the width of the normal distributions used for proposing new values of beta and gamma
# Smaller "sd_proposals" mean more steps accepted but slower movement through parameter space
# Larger "sd_proposals" mean lots of rejection but broader exploration (hopping around a lot)

# Initialisation:
bchain = rep(0, nsamples); # Vector column of "nsample" zeros - values of beta
gchain = rep(0, nsamples); # Vector column of "nsample" zeros - values of gamma
LLchain = rep(0, nsamples); # Log-likelihood

pc = c(5, 1); # Current (c) value of parameters [beta,gamma]
# Initial values should be computed from the prior, but it doesn't matter
# as they don't influence the results if the chain has been run long enough
logpriorg_c = (shape-1)*log(pc[2]) - pc[2]/scale;  
LLc = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pc); # Current value of the log-likelihood (LL), i.e. calculated in current parameter value

accept_count = 0; # This is to count how many times I accept if I propose both beta and gamma together (block proposal)
accept_count_b = 0; # This is to count how many times I accept proposed values of beta
accept_count_g = 0; # This is to count how many times I accept proposed values of gamma

for (s in 1:nsamples){
    bchain[s] = pc[1]; # Store current values of beta
    gchain[s] = pc[2]; # Store current values of gamma
    LLchain[s] = LLc; # Store current value of the log-likelihood
    if (block_proposal){ # If I do block updating (both beta and gamma at the same time
        for (t in 1:(2*thinning)){ # This is to make the same number of moves in block and separate updates
            pp = mvrnorm(1,pc,diag(sd_proposals)); # Propose both parameters from a multivariate normal distribution
#             pp = mvrnorm(1,pc,diag(sd_proposals)); # This is a lot more efficient
            if (all(pp>0)) { # Continue only if both parameters are > 0, else reject immediately
                logpriorg_p = (shape-1)*log(pp[2]) - pp[2]/scale; # This is the log-prior for the new value of gamma (up to a constant)
#                 logpriorg_p = log(gampdf(pp[2],shape,scale)); # Just to check I've done th maths right!
                LLp = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pp); # Calculate log-likelihood of proposal
                logalpha = LLp + logpriorg_p - LLc - logpriorg_c; # log of probability of acceptance
                if (logalpha > log(runif(1))) { # Accept
                    pc = pp; # Save proposed value as new value
                    logpriorg_c = logpriorg_p; # Save value of log-prior of newly accepted gamma
                    LLc = LLp; # Save log-likelihood of newly accepted value
                    accept_count = accept_count + 1; # Record the acceptance event
                } # else leave pc and LLc as they were before
            }
        }
    }
    else { #If I don't do block updating, but I update one parameter at a time
        for (t in 1:thinning) {
            # Update beta
            pb = rnorm(1, pc[1],sd_proposals[1]); # Propose a value of beta from (1-dim) normal 
            if (pb > 0) { # Otherwise, reject immediately
                pp = c(pb, pc[2]); # save new value in a vector of parameters (old gamma stays the same)
                # gamma hasn't changed, so no need to compute the new value for the prior
                LLp = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pp); # New log-likelihood
                logalpha = LLp - LLc; # log acceptance probability (only beta changed, and beta has flat prior)
                if (logalpha > log(runif(1))) { # Accept
                    pc = pp; # Save proposed value as new value
                    LLc = LLp; # Save log-likelihood of newly accepted value
                    accept_count_b = accept_count_b + 1; # Record the acceptance event
                } # else leave pc and LLC as they were before
            }
            # Update gamma
            pg = rnorm(1, pc[2],sd_proposals[2]);
            if (pg > 0) { # Otherwise, reject immediately
                pp = c(pc[1], pg); # vector of proposed parameters (beta hasn't changed)
                logpriorg_p = (shape-1)*log(pp[2]) - pp[2]/scale;  # new value for the prior
                LLp = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pp);
                logalpha = LLp + logpriorg_p - LLc - logpriorg_c; # log-prob of acceptance (this time I need the prior for gamma)
                if (logalpha > log(runif(1))) { # Accept
                    pc = pp; # Save proposed value as new value
                    logpriorg_c = logpriorg_p; # Save value of log-prior of newly accepted gamma
                    LLc = LLp; # Save log-likelihood of newly accepted value
                    accept_count_g = accept_count_g + 1; # Record the acceptance event
                } # else leave pc and LLC as they were before
            }
        }
    }        
    if (s %% 500 == 0) cat(paste('Loop ', s, ' of ', nsamples, ' completed. \n'))
}

## PLOTTING ######

# To run only this section in R, click  anywhere below the "##" just
# above here and press "Code -> Run Region -> Run Code Section" (or, on Mac, alt+cmd+T)

# Parameters for post-MCMC analysis of the data
burnin = 1; # Define length of burnin period, i.e. how many samples to discard
nbins = 100; # Numer of bins for the histogram

# Trace plots
par(mfrow=c(3,1)) # Divide figure in 3 parts (3 rows, 1 column)
plot(burnin:nsamples, LLchain[burnin:nsamples], typ="l", col = "black",
     xlab = "samples", ylab = "LL") # Plot log-likelihood (excluding burnin) against sample number
plot(burnin:nsamples,bchain[burnin:nsamples], typ="l", col = "red", 
     xlab = "Samples", ylab = "beta") # Plot sampled values of beta (excluding burnin)
plot(burnin:nsamples,gchain[burnin:nsamples], typ="l", col = "blue",
     xlab = 'Samples', ylab = 'gamma') # Plot sampled values of gamma (excluding burnin)

# Autocorrelation
par(mfrow=c(2,1)) # Divide figure in 2 parts (2 rows, 1 column)
acf(bchain[burnin:nsamples],min(nsamples-1,100)); # Draw the autocorrelation plot for beta
acf(gchain[burnin:nsamples],min(nsamples-1,100)); # Draw the autocorrelation plot for gamma
# This requires the function acf (built-in to R)

# Posterior distribution
if (block_proposal) { # If I do plock proposal, I don't count separate acceptance rates, 
    # but the pictures aske for them in the title, so I make them both equal to the common value
    accept_count_b = accept_count;
    accept_count_g = accept_count;
}

par(mfrow=c(2,1)) # Divide figure in 2 parts (2 rows, 1 column)
h = hist(bchain[burnin:nsamples],nbins,freq=FALSE, main=NULL); # Plot the histogram of p 
# (excluding burnin), with "nbins" bins, and normalise it (rather than counting total numbers in each bin)
# The prior for beta is improper (i.e., because beta is defined in an
# infinitely long set, i.e. [0,+infinity), it is 0 everywhere), so I don't plot it
# Hard to compute the posterior analytically, so don't plot posterior
ux <- unique(bchain[burnin:nsamples])
md = ux[which.max(tabulate(match(bchain[burnin:nsamples], ux)))] # calculate mode
title(paste('Mean = ', round(1e5*mean(bchain[burnin:nsamples]))/1e5, 
            '; mode = ', round(1e5*md)/1e5,
            '; sd = ', round(1e5*sd(bchain[burnin:nsamples]))/1e5,
            '; acceptance rate = ', round(1e4*100*accept_count_b/niters)/1e4, '%'))
# Write title, using "num2str" to transform numbers into strings, "mean", "mode" 
# and "std" to calculate mean, mode and standard deviation of the sample of beta,
# and calculating the acceptance probability
# "..." means: continue command on next line

hist(gchain[burnin:nsamples],nbins,freq=FALSE, main=NULL); # Plot the histogram of p 
# (excluding burnin), with "nbins" bins, and normalise it (rather than counting total numbers in each bin)
# If you want to zoom the x-axis only where the histogram is, uncomment the
# following 2 lines and comment the third one
xl = range(gchain); # Get the limits for x-axis set automatically by matlab, to plot the prior
x = seq(xl[1], xl[2], 0.01); # Create a vector of value for the x axis, with step 0.01
lines(x,dgamma(x, shape, 1/scale), typ="l", col = "magenta"); # Plot the prior for gamma (not sure why it has to be 1/scale)
# Hard to compute the posterior analytically, so don't plot
ux <- unique(gchain[burnin:nsamples])
md = ux[which.max(tabulate(match(gchain[burnin:nsamples], ux)))] # calculate mode
title(paste('Mean = ', round(1e5*mean(gchain[burnin:nsamples]))/1e5, 
            '; mode = ', round(1e5*md)/1e5,
            '; sd = ', round(1e5*sd(gchain[burnin:nsamples]))/1e5,
            '; acceptance rate = ', round(1e4*100*accept_count_g/niters)/1e4, '%'))
# Write title, using "num2str" to transform numbers into strings, "mean", "mode" 
# and "std" to calculate mean, mode and standard deviation of the sample of gamma,
# and calculating the acceptance probability
# "..." means: continue command on next line

par(mfrow=c(1,1)) # Divide figure in 3 parts (3 rows, 1 column)
scatterhist(bchain[burnin:nsamples], gchain[burnin:nsamples],
     xlab = "beta", ylab = "gamma") # Plot correlation between two parameters

toc = proc.time() - tic # finish timer
toc # show time

