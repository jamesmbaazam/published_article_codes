Gamma_meansd2shapescale = function(m,s){
  
# This function changes the parameterisation of the Gamma distribution from
# mean and standard devation to the shape and scale parameters used by matlab
#
# In R, 
# a is the shape parameter (alpha) and 
# b is the scale parameter, defined as b = 1/lambda
#
# mean = alpha/lambda = a * b
# var = alpha/lambda^2 = a * b^2
# var = mean * b --> b = var / mean = s^2 / m
# a = mean / b = mean^2 / var = m^2 / s^2

shape = m^2/s^2;
scale = s^2/m;

vals = list(shape, scale)
return(vals)

}