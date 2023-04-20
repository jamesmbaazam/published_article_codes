# This script plots the histogram for the binomial distribution

# Parameters (to play with)
n = 10; # number of coin tosses
p = 0.5; # probability of head

# Initialisation
x = 0:n; # List of possible values for the number of successes (row vector, x-axis)

# Main program
y = dbinom(x,n,p);
# Matlab allows x to be a vector, and spits out y of the same dimensions as x.
# Otherwise a for loop would be needed
y # Display y

# Plot
barplot(y, xlab='Number of successes', ylab='Probability')
title(main = 'Binomial distribution')