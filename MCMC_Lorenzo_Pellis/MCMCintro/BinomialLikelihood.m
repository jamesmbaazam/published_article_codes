% This script plots the histogram for the many binomial distributions at
% the same time

% Parameters (to play with)
n = 10; % number of coin tosses
xvec = [ 0 1 3 5 8 ]; % vector of x values observed
lx = length(xvec); % number of different values of x = number of distributions to plot
dp = 0.001; % step in the values of p
pvec = (0:dp:1)'; % vector of values of p (in column)
lp = length(pvec);

% Initialisation
l = zeros(lp,lx); % one column per distribution

% Main program
for ip = 1:lp % index for p ranges from 1 to lp
    p = pvec(ip);
    l(ip,:) = binopdf(xvec,n,p); % Store results in the ip's row
    % Matlab allows x to be a vector, and spits out y of the same dimensions as x.
    % Otherwise a for loop would be needed
end

% Plot
figure(2)
plot(pvec,l,'Linewidth',2)
xlabel('Probability of success p')
ylabel('Likelihood')
title('Binomial Likelihood (10 trials)')
leg = cell(1,lx);
for ix = 1:lx
    leg{ix} = [num2str(xvec(ix)),' successes'];
end
legend(leg)