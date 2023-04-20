% MCMC code to estimate transmission and recovery rates from the final size
% of an SIR epidemic in a large population. The final size is measured by
% sampling "n_tested" individuals (with a perfect test) and detecting that 
% "n_pos" have been infected

close all; % Close all figures left open from previous sessions
clearvars; % Clear all variables to start fresh!

% Data given
n_tested = 100;
n_pos = 80;

% Propose new parameters in block or one at a time?
block_proposal = false; 
% if true, propose both parameters from a bivariate normal distribution
% if false, propose one parameter at a time

% Priors:
% For beta: flat, improper prior on [0,+infinity)
% For gamma: gamma-distributed with mean m_hy_gamma days and sd sd_hy_gamma days
m_hy_gamma = 4; % hyper-parameter for shape of Gamma distribution for gamma
sd_hy_gamma = 1; % hyper-parameter for scale of Gamma distribution of gamma
% Matlab parameterises the gamma distribution in terms of shape and scale parameters
[ shape, scale ] = Gamma_meansd2shapescale(m_hy_gamma,sd_hy_gamma); % Find corresponding shape and scale

% MCMC parameters
niters = 10000; % Number of iterations
thinning = 1; % Save only 1 iteration every "thinning"
nsamples = floor(niters/thinning); % Effective number of samples that will be saved
% "floor" truncates to integer immediately below, i.e. floor(3.7623) = 3
sd_proposals = [ 2, 1 ]; % This is the width of the normal distributions used for proposing new values of beta and gamma
% Smaller "sd_proposals" mean more steps accepted but slower movement through parameter space
% Larger "sd_proposals" mean lots of rejection but broader exploration (hopping around a lot)

% Initialisation:
bchain = zeros(nsamples,1); % Vector column of "nsample" zeros - values of beta
gchain = zeros(nsamples,1); % Vector column of "nsample" zeros - values of gamma
LLchain = zeros(nsamples,1); % Log-likelihood

pc = [ 5, 1 ]; % Current (c) value of parameters [beta,gamma]
% Initial values should be computed from the prior, but it doesn't matter
% as they don't influence the results if the chain has been run long enough
logpriorg_c = (shape-1)*log(pc(2)) - pc(2)/scale;  
LLc = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pc); % Current value of the log-likelihood (LL), i.e. calculated in current parameter value

accept_count = 0; % This is to count how many times I accept if I propose both beta and gamma together (block proposal)
accept_count_b = 0; % This is to count how many times I accept proposed values of beta
accept_count_g = 0; % This is to count how many times I accept proposed values of gamma

for s = 1:nsamples
    bchain(s) = pc(1); % Store current values of beta
    gchain(s) = pc(2); % Store current values of gamma
    LLchain(s) = LLc; % Store current value of the log-likelihood
    if block_proposal % If I do block updating (both beta and gamma at the same time
        for t = 1:(2*thinning) % This is to make the same number of moves in block and separate updates
            pp = mvnrnd(pc,diag(sd_proposals)); % Propose both parameters from a multivariate normal distribution
%             pp = mvnrnd(pc,diag(sd_proposals)); % This is a lot more efficient
            if all(pp>0) % Continue only if both parameters are > 0, else reject immediately
                logpriorg_p = (shape-1)*log(pp(2)) - pp(2)/scale; % This is the log-prior for the new value of gamma (up to a constant)
%                 logpriorg_p = log(gampdf(pp(2),shape,scale)); % Just to check I've done th maths right!
                LLp = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pp); % Calculate log-likelihood of proposal
                logalpha = LLp + logpriorg_p - LLc - logpriorg_c; % log of probability of acceptance
                if logalpha > log(rand) % Accept
                    pc = pp; % Save proposed value as new value
                    logpriorg_c = logpriorg_p; % Save value of log-prior of newly accepted gamma
                    LLc = LLp; % Save log-likelihood of newly accepted value
                    accept_count = accept_count + 1; % Record the acceptance event
                end % else leave pc and LLc as they were before
            end
        end
    else % If I don't do block updating, but I update one parameter at a time
        for t = 1:thinning
            % Update beta
            pb = normrnd(pc(1),sd_proposals(1)); % Propose a value of beta from (1-dim) normal 
            if pb > 0 % Otherwise, reject immediately
                pp = [ pb, pc(2) ]; % save new value in a vector of parameters (old gamma stays the same)
                % gamma hasn't changed, so no need to compute the new value for the prior
                LLp = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pp); % New log-likelihood
                logalpha = LLp - LLc; % log acceptance probability (only beta changed, and beta has flat prior)
                if logalpha > log(rand) % Accept
                    pc = pp; % Save proposed value as new value
                    LLc = LLp; % Save log-likelihood of newly accepted value
                    accept_count_b = accept_count_b + 1; % Record the acceptance event
                end % else leave pc and LLC as they were before
            end
            % Update gamma
            pg = normrnd(pc(2),sd_proposals(2));
            if pg > 0 % Otherwise, reject immediately
                pp = [ pc(1), pg ]; % vector of proposed parameters (beta hasn't changed)
                logpriorg_p = (shape-1)*log(pp(2)) - pp(2)/scale;  % new value for the prior
                LLp = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pp);
                logalpha = LLp + logpriorg_p - LLc - logpriorg_c; % log-prob of acceptance (this time I need the prior for gamma)
                if logalpha > log(rand) % Accept
                    pc = pp; % Save proposed value as new value
                    logpriorg_c = logpriorg_p; % Save value of log-prior of newly accepted gamma
                    LLc = LLp; % Save log-likelihood of newly accepted value
                    accept_count_g = accept_count_g + 1; % Record the acceptance event
                end % else leave pc and LLC as they were before
            end
        end
    end        
    disp(['Loop ',num2str(s),' of ',num2str(nsamples),' completed!']); % Display completed loop
end

%%

% To run only this section in matlab, click  anywhere below the "%%" just
% above here (this section's background turns yellow) and press "Run Section"

% Parameters for post-MCMC analysis of the data
burnin = 1; % Define length of burnin period, i.e. how many samples to discard
nbins = 100; % Numer of bins for the histogram

% Trace plots
scrsz = get(groot,'ScreenSize'); % get screensize for shaping figures
figure(1)
set(gcf,'Position',[ scrsz(3)/2 scrsz(4)*3/4 scrsz(3)/2 scrsz(4)*3/4 ]) % Create new figure with appropriate position and size
subplot(3,1,1) % Divide figure in 3 parts (3 rows, 1 column) and select part 1 (top)
plot(burnin:nsamples,LLchain(burnin:nsamples),'k') % Plot log-likelihood (excluding burnin) against sample number
xlabel('Samples') % write label to x-axis
ylabel('LL') % write label to y-axis

subplot(3,1,2) % Focus on part 2 of the figure (middle)
plot(burnin:nsamples,bchain(burnin:nsamples),'r') % Plot sampled values of beta (excluding burnin)
xlabel('Samples')
ylabel('beta')

subplot(3,1,3) % Focus on part 3 of the figure (bottom)
plot(burnin:nsamples,gchain(burnin:nsamples),'b') % Plot sampled values of gamma (excluding burnin)
xlabel('Samples')
ylabel('gamma')

% Autocorrelation
figure(2)
set(gcf,'Position',[ 1 scrsz(4)/10 scrsz(3)/3 scrsz(4)/2 ]) % Create new figure with appropriate position and size
subplot(2,1,1) % top subplot
acf(bchain(burnin:nsamples),min(nsamples-1,100)); % Draw the autocorrelation plot for beta
subplot(2,1,2) % bottom subplot
acf(gchain(burnin:nsamples),min(nsamples-1,100)); % Draw the autocorrelation plot for gamma
% This requires the function acf (downloaded from the web)

% Postrior distribution
if block_proposal % If I do plock proposal, I don't count separate acceptance rates, 
    % but the pictures aske for them in the title, so I make them both equal to the common value
    accept_count_b = accept_count;
    accept_count_g = accept_count;
end
figure(3) % Create a new figure
set(gcf,'Position',[ 1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2 ]) % Create new figure with deisred position and size
clf % Clean figure, i.e. cancel anything that was left on this figure and start fresh
subplot(2,1,1) % top subplot
h = histogram(bchain(burnin:nsamples),nbins,'Normalization','pdf'); % Plot the histogram of p 
% (excluding burnin), with "nbins" bins, and normalise it (rather than counting total numbers in each bin)
% The prior for beta is improper (i.e., because beta is defined in an
% infinitely long set, i.e. [0,+infinity), it is 0 everywhere), so I don't plot it
% Hard to compute the posterior analytically, so don't plot posterior
[ ~, index ] = max(h.Values); % this is just a hack to compute the mode of the histogram
title(['Mean = ',num2str(mean(bchain(burnin:nsamples)),'%0.3g'),'; mode = ',num2str(h.BinEdges(index)+h.BinWidth*0.5,'%0.3g'),...
    '; sd = ',num2str(std(bchain(burnin:nsamples)),'%0.3g'),'; acceptance rate = ',num2str(accept_count_b/niters,'%.2f')]);
% Write title, using "num2str" to transform numbers into strings, "mean", "mode" 
% and "std" to calculate mean, mode and standard deviation of the sample of beta,
% and calculating the acceptance probability
% "..." means: continue command on next line

subplot(2,1,2)
h = histogram(gchain(burnin:nsamples),nbins,'Normalization','pdf'); % Plot the histogram of p 
% (excluding burnin), with "nbins" bins, and normalise it (rather than counting total numbers in each bin)
% If you want to zoom the x-axis only where the histogram is, uncomment the
% following 2 lines and comment the third one
xl = xlim; % Get the limits for x-axis set automatically by matlab, to plot the prior
x = xl(1):0.01:xl(2); % Create a vector of value for the x axis, with step 0.01
hold on; % Keep what is in figure when plotting other things (following lines)
plot(x,gampdf(x,shape,scale),'m'); % Plot the prior for gamma
% Hard to compute the posterior analytically, so don't plot
[ ~, index ] = max(h.Values); % hack to compute the mode of the histogram
title(['Mean = ',num2str(mean(gchain(burnin:nsamples)),'%0.3g'),'; mode = ',num2str(h.BinEdges(index)+h.BinWidth*0.5,'%0.3g'),...
    '; sd = ',num2str(std(gchain(burnin:nsamples)),'%0.3g'),'; acceptance rate = ',num2str(accept_count_g/niters,'%.2f')]);
% Write title, using "num2str" to transform numbers into strings, "mean", "mode" 
% and "std" to calculate mean, mode and standard deviation of the sample of gamma,
% and calculating the acceptance probability
% "..." means: continue command on next line

figure(4)
set(gcf,'Position',[ scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/1.6 ]) % Create new figure
plotmatrix([bchain(burnin:nsamples),gchain(burnin:nsamples)]);
% plotmatrix is a useful command for the fancy graph to explore correlations between the parameters

