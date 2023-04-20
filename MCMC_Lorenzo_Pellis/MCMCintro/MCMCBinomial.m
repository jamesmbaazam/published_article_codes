% MCMC code to estimate the probability of a coin giving head

close all; % Close all figures left open from previous sessions
clearvars; % Clear all variables to start fresh!

% Data given
n_tosses = 100;
n_heads = 79;

% If you want to generate other data, uncomment the following:
% p_true = 0.8;
% n_heads = binornd(n_tosses,p_true);

% Prior (uniform) - don't need to specify anything...

% MCMC parameters
niters = 10000; % Number of iterations
thinning = 1; % Save only 1 iteration every "thinning"
nsamples = floor(niters/thinning); % Effective number of samples that will be saved
% "floor" truncates to integer immediately below, i.e. floor(3.7623) = 3
sd_proposal = 0.5; % This is the width of the normal distribution used for the proposal
% Smaller "sd_proposal" mean more steps accepted but slower movement through parameter space
% Larger "sd_proposal" mean lots of rejection but broader exploration (hopping around a lot)

% Initialisation:
pchain = zeros(nsamples,1); % Vector column of "nsample" zeros
LLchain = zeros(nsamples,1); % Ditto

pc = 0.2; % Current (c) value of parameter p, drawn from the prior 
% I can also fix it as I want, given it bears no influence on the chain if
% it has "converged" (run long enough), as it should have...
LLc = log(binopdf(n_heads,n_tosses,pc)); % Current value of the log-likelihood (LL), i.e. calculated in current parameter value
% because the binomial coefficient is constant (i.e. indpendent of pc), we could also directly calculate this analytically as:
% LLc = n_heads * log(pc) + (n_tosses-n_heads)*log(1-pc)

accept_count = 0; % This is to count how many times I accept

for s = 1:nsamples
    pchain(s) = pc; % Store current values of parameter
    LLchain(s) = LLc; % Store current value of the log-likelihood
    for t = 1:thinning
        pp = normrnd(pc,sd_proposal); % Proposed (p) new value for the parameter
        if ( pp > 0 ) && ( pp < 1 )
            LLp = log(binopdf(n_heads,n_tosses,pp)); % Compute the new log-likelihood
            logalpha = LLp - LLc; % Compute the log  of the probability of acceptance
            % logalpha has a particularly simple form because: 
            % 1) The prior is uniform (values cancel out), and 
            % 2) the proposal distribution is symmetric (values cancel out)
            if logalpha >= 0 % i.e. if alpha >= 1, accept for sure
                pc = pp; % Save proposed value as new value
                LLc = LLp; % Save log-likelihood of newly accepted value
                accept_count = accept_count + 1; % Record the acceptance event
            else % if alpha < 1
                r = rand; % Generate a random number to decide whether to accept the new value or not
                if log(r) < logalpha % If r < alpha, accept (if in doubt, < is correct, because if alpha = 1 I would always accept)
                    pc = pp; % Save proposed value as new value
                    LLc = LLp; % Save log-likelihood of newly accepted value
                    accept_count = accept_count + 1; % Record the acceptance event
                end
            end
        end
    end
    disp(['Loop ',num2str(s),' of ',num2str(nsamples),' completed!']); % Display completed loop
end

%%

% To run only this section in matlab, click  anywhere below the "%%" just
% above here (this section's background turns yellow) and press "Run Section"

% Parameters for post-MCMC analysis of the data
burnin = 100; % Define length of burnin period, i.e. how many samples to discard
nbins = 100; % Numer of bins for the histogram

% Trace plots
figure(1) % Create new figure
subplot(2,1,1) % Divide figure in 2 parts (2 rows, 1 column) and select part 1 (top)
plot(burnin:nsamples,LLchain(burnin:nsamples),'k') % Plot log-likelihood (excluding burnin) against sample number
xlabel('Samples') % write label to x-axis
ylabel('LL') % write label to y-axis

subplot(2,1,2) % Focus on part 2 of the figure (bottom)
plot(burnin:nsamples,pchain(burnin:nsamples),'b') % Plot sampled values of p (excluding burnin)
xlabel('Samples')
ylabel('p')

% Autocorrelation
figure(2) % Create a new figure
acf(pchain(burnin:nsamples),min(nsamples-1,100)); % Draw the autocorrelation plot
% This requires the function acf (downloaded from the web)

% Postrior distribution
figure(3) % Create a new figure
clf % Clean figure, i.e. cancel anything that was left on this figure and start fresh
histogram(pchain(burnin:nsamples),nbins,'Normalization','pdf') % Plot the histogram of p 
% (excluding burnin), with "nbins" bins, and normalise it (rather than counting total numbers in each bin)
% If you want to zoom the x-axis only where the histogram is, uncomment the
% following 2 lines and comment the third one
% xl = xlim; % Get the limits for x-axis set automatically by matlab 
% x = xl(1):0.01:xl(2); % Create a vector of value for the x axis, with step 0.01
x = 0:0.01:1; % Create a vector of value for the x axis from 0 to 1, with step 0.01
hold on; % Keep what is in figure when plotting other things (following lines)
plot(x,ones(size(x)),'r'); % Plot the prior, i.e. horizontal line at 1
plot(x,betapdf(x,1+n_heads,1+n_tosses-n_heads),'k'); % Plot the exact value of the posterior (which we know only because this is a simple case!)
title(['Mean = ',num2str(mean(pchain(burnin:nsamples))),'; mode = ',num2str(mode(pchain(burnin:nsamples))),...
    '; sd = ',num2str(std(pchain(burnin:nsamples))),'; acceptance rate = ',num2str(accept_count/niters,'%.2f')]);
% Write title, using "num2str" to transform numbers into strings, "mean", "mode" 
% and "std" to calculate mean, mode and standard deviation of the sample of p,
% and calculating the acceptance probability (written to 2 decimal places)
% "..." means: continue command on next line



