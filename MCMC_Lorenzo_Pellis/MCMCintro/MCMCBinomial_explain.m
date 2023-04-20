% MCMC code to estimate the probability p of a coin giving heads

close all;
clearvars;

height = 3;
minx = 0;
maxx = 1;
dx = 0.001;
x = (minx:dx:maxx)';
lx = length(x);
minxw = -0.2;
maxxw = 1.2;
xw = (minxw:dx:maxxw)';
lxw = length(xw);
rng(13);
scrsz = get(groot,'ScreenSize');
figure('Position',[ scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2 ])
plot([minx minx],[0 height],':k')
hold on
plot([maxx maxx],[0 height],':k')
xlim([minxw maxxw])
ylim([0 height])
set(gca,'XTick',minxw:0.1:maxxw);
pause;

% Generate the data
p_true = 0.8;
n_tosses = 10;
n_heads = binornd(n_tosses,p_true);

% Explore the data
p_hat = n_heads / n_tosses;

% Prior distribution - uniform in [0,1) (and 0 outside)
y_prior = zeros(size(xw));
y_prior((xw>=0)&(xw<=1)) = 1;
hb = plot(xw,y_prior,'r');
hbt = text(maxx,1.1,'Prior distribution','Color','r');
hbtA = text(0.1,0.5,'Area = 1','color','r');
pause;
set(hbtA,'Visible','off')

% Likelihood
y_likelihood = zeros(size(x));
for ix = 1:length(x)
    y_likelihood(ix) = binopdf(n_heads,n_tosses,x(ix));
end
hl = plot(x,y_likelihood,'g');
hlt = text(p_hat+0.1,max(y_likelihood),'Likelihood','Color','g','VerticalAlignment','bottom');
area = trapz(x,y_likelihood);
hltA = text(p_hat,max(y_likelihood)/2,['Area = ',num2str(area)],'Color','g','HorizontalAlignment','center');
pause;
y_posterior = y_likelihood / area;
ha = plot(x,y_posterior,'k');
hat = text(p_hat+0.1,max(y_posterior),'Posterior','Color','k');
hatA = text(p_hat,max(y_posterior)/2,['Area = 1'],'Color','k','HorizontalAlignment','center');
pause;
set(hl,'Visible','off')
set(hlt,'Visible','off')
set(hltA,'Visible','off')
set(hatA,'Visible','off')
% y_posterior_analytical = betapdf(x,1+n_heads,1+n_tosses-n_heads);
% haa = plot(x,y_posterior_analytical,'--g');

% MCMC parameters
niters = 100;
sigma_proposal = 0.1;

pchain = zeros(niters,1);
Lchain = zeros(niters,1);
hpchain = zeros(niters,1);

harr = 0.2;
pc = rand; % starting value of p chosen from the prior
Lc = binopdf(n_heads,n_tosses,pc);
hc = plot(pc,0,'vm');
hcarr = plot([pc pc], [0 harr],'m');
pause;
htc = text(pc,2*harr,{['p_c = ',num2str(pc)];['L_c = ',num2str(Lc)]},'color','m','HorizontalAlignment','center');
pause;
acc_count = 0;

for it = 1:niters
    pchain(it) = pc;
    Lchain(it) = Lc;
    
    % Draw proposal distribution
    xp = (ceil((pc-3*sigma_proposal)/dx)*dx):dx:(floor((pc+3*sigma_proposal)/dx)*dx);
    yp = normpdf(xp,pc,sigma_proposal)/5;
    hpd = plot(xp,yp,'m');
    hpdt = text(pc+0.1,max(yp),{'Proposal distribution';'(not to scale)'},'color','m');
    pause;
    pp = normrnd(pc,sigma_proposal);
    hpp = plot(pp,0,'vb');
    hparr = plot([pp pp], [0 harr],'--b');
    pause;
    htp = text(0,2,{['p_p = ',num2str(pp)]},'color','b');
    pause;
    if ( pp < 0 ) || ( pp >= 1 )
        % Don't do anything
        hte = text(0,1.7,{'The proposed value is outside the prior range';'Reject, i.e. keep old value'},'color','b');
    else
        Lp = binopdf(n_heads,n_tosses,pp);
        alpha = Lp / Lc;
        if alpha > 1
            hte = text(0,1.5,{['L_p = ',num2str(Lp)];['alpha = L_p/L_c = ',num2str(alpha)];'alpha > 1 --> accept'},'color','b');
            pc = pp;
            Lc = Lp;
        else
            r = rand;
            if r < alpha
                hte = text(0,1.5,{['L_p = ',num2str(Lp)];['alpha = L_p/L_c = ',num2str(alpha)];['Generate random number r = ',num2str(r)];'r < alpha --> accept'},'color','b');
                pc = pp;
                Lc = Lp;
            else
                % Don't do anything
                hte = text(0,1.5,{['L_p = ',num2str(Lp)];['alpha = L_p/L_c = ',num2str(alpha)];['Generate random number r = ',num2str(r)];'r > alpha --> reject'},'color','b');
            end
        end
    end
    pause;
    set(hc,'Visible','off')
    set(hcarr,'Visible','off')
    set(htc,'Visible','off')
    set(hpd,'Visible','off')
    set(hpdt,'Visible','off')
    set(hparr,'Visible','off')
    set(hpp,'Visible','off')
    set(htp,'Visible','off')
    set(hte,'Visible','off')
    hpchain(it+1) = plot(pc,0,'*c');
    pause;

    disp(['Loop ',num2str(it),' of ',num2str(niters),' completed!']);

    hc = plot(pc,0,'vm');
    hcarr = plot([pc pc], [0 harr],'m');
    pause;
    htc = text(pc,2*harr,{['p_c = ',num2str(pc)];['L_c = ',num2str(Lc)]},'color','m','HorizontalAlignment','center');

end

