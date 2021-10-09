% This script plots the histogram for the many binomial distributions at
% the same time

% Parameters (to play with)
n = 10; % number of coin tosses
pvec = [ 0.1 0.5 0.7 ]; % vector for probabilities of head
lp = length(pvec); % number of distributions plotted simultaneously

% Initialisation
x = 0:n; % List of possible values for the number of successes (row vector, x-axis)
y = zeros(lp,n+1); % one row per distribution, each of them with n+1 possible number of heads

% Main program
for ip = 1:lp % index for p ranges from 1 to lp
    p = pvec(ip);
    y(ip,:) = binopdf(x,n,p); % Store results in the ip's row
    % Matlab allows x to be a vector, and spits out y of the same dimensions as x.
    % Otherwise a for loop would be needed
end
y % Display y

% Plot
figure(1)
bar(x,y') % matlab prefers if each group is in one column, rather than 1 row, so I transpose y
xlabel('Number of successes')
ylabel('Probability')
title('Binomial distributions out of 10 trials')
leg = cell(1,lp);
for ip = 1:lp
    leg{ip} = ['p = ',num2str(pvec(ip))];
end
legend(leg)
