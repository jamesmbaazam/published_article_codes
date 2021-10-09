LL_EpidemicFinalSizeLargePop = function(n_pos,n_tested,pc){
  
  require(utils)
  # pc = [beta,gamma]
  b = pc[1]
  g = pc[2]
  
  R0 = b/g;
  if (R0 <= 1) z = 0
  else{
    
    f = function (s, R0) 1 - s - exp(-R0*s)
    y = uniroot(f, c(1e-8, 1), tol = 0.0001, R0 = b/g) # get value of s for f(s) = 0
    z = y$root
    if (is.nan(z)){
      beep;
      disp('Likelihood spit out NaN - problems might occur')
      pause;
    }
  }
  
  LL = n_pos * log(z) + (n_tested-n_pos) * log(1-z); # likelihood
  return(LL) # output LL from function
  
}