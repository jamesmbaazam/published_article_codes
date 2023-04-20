function LL = LL_EpidemicFinalSizeLargePop(n_pos,n_tested,pc)

% pc = [beta,gamma]
b = pc(1);
g = pc(2);

R0 = b/g;
if R0 <= 1
    z = 0;
else
    funz = @(s) ( 1 - s - exp(-R0*s) );
    z = fzero(funz,[1e-8,1]);
    if isnan(z)
        beep;
        disp('Likelihood spit out NaN - problems might occur')
        pause;
    end
end

LL = n_pos * log(z) + (n_tested-n_pos) * log(1-z);