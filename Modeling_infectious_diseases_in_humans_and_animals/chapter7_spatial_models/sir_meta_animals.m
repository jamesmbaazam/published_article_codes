From: <Saved by Blink>
Snapshot-Content-Location: http://www.modelinginfectiousdiseases.org/
Subject: Modeling Infectious Diseases in Humans and Animals, by Keeling & Rohani
Date: Wed, 12 Feb 2020 14:33:56 -0000
MIME-Version: 1.0
Content-Type: multipart/related;
	type="text/html";
	boundary="----MultipartBoundary--O7eOMlZSlJC40B9HEldG57t80vTSAZsSRTgVi59Jih----"


------MultipartBoundary--O7eOMlZSlJC40B9HEldG57t80vTSAZsSRTgVi59Jih----
Content-Type: text/html
Content-ID: <frame-902DFE8F8E4AA2B90B2BA07557E43BDC@mhtml.blink>
Content-Transfer-Encoding: quoted-printable
Content-Location: http://www.modelinginfectiousdiseases.org/

<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Frameset//EN" "http://www.w3.o=
rg/TR/html4/frameset.dtd"><html lang=3D"en"><head><meta http-equiv=3D"Conte=
nt-Type" content=3D"text/html; charset=3DUTF-8"><title>Modeling Infectious =
Diseases in Humans and Animals, by Keeling &amp; Rohani</title>
<meta name=3D"Keywords" content=3D"Modeling Infectious Diseases in humans a=
nd animals, Princeton University Press,"><meta name=3D"Author" content=3D"M=
att Keeling &amp; Pejman Rohani."><meta name=3D"Description" content=3D"The=
se web pages contain all the programs labelled in the book &quot;Modeling I=
nfectious Diseases in Humans and Animals&quot;. They are generally availabl=
e as C++, Fortran and Matlab files."><link rel=3D"shortcut icon" href=3D"ht=
tp://go.warwick.ac.uk/ModelingInfectiousDiseases/favicon.ico"></head>
<frameset rows=3D"100%" data-gr-c-s-loaded=3D"true">
<frame title=3D"http://go.warwick.ac.uk/ModelingInfectiousDiseases" src=3D"=
cid:frame-A68C0874F133FB69827CA7B11B8F8D46@mhtml.blink" name=3D"mainframe" =
frameborder=3D"0" noresize=3D"noresize" scrolling=3D"auto">
<noframes>Sorry, you don"t appear to have frame support.
Go here instead - <a href=3D"http://go.warwick.ac.uk/ModelingInfectiousDise=
ases">Modeling Infectious Diseases in Humans and Animals, by Keeling & Roha=
ni</a></noframes>
</frameset></html>
------MultipartBoundary--O7eOMlZSlJC40B9HEldG57t80vTSAZsSRTgVi59Jih----
Content-Type: text/html
Content-ID: <frame-A68C0874F133FB69827CA7B11B8F8D46@mhtml.blink>
Content-Transfer-Encoding: quoted-printable
Content-Location: http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter7/Program_7.1/Program_7_1.m

<html><head><meta http-equiv=3D"Content-Type" content=3D"text/html; charset=
=3Dwindows-1252"></head><body><pre style=3D"word-wrap: break-word; white-sp=
ace: pre-wrap;">function [t,X,Y] =3D Program_7_1(n,beta,gamma,nu,mu,m,X0,Y0=
,MaxTime)
%
%=20
%
% Program_7_1( n, beta, gamma, nu, mu, m, X0, Y0, MaxTime)
%      This is the MATLAB version of program 7.1 from page 241 of=20
% "Modeling Infectious Disease in humans and animals"=20
% by Keeling &amp; Rohani.
%=20
% It is the SIR epidemic in a metapopulation with "animal-like" movements
% of infected or susceptible individuals across the network.
%=20
% The default model employs global coupling between ALL subpopulations; an
% alterative linear arrangement, with nearest-neighbour coupling can be
% achieved by:
%=20
% Program_7_1(n,beta,gamma,nu,mu,m*(diag(ones(1,n-1),1)+diag(ones(1,n-1),-1=
)),X0,Y0,MaxTime)
%

% Sets up default parameters if necessary.
if nargin =3D=3D 0
   n=3D5;
   beta=3D1.0*ones(n,1);
   gamma=3D0.1*ones(n,1);
   nu=3D0.0001*ones(n,1);
   mu=3D0.0001*ones(n,1);
   X0=3D0.1*ones(n,1);
   Y0=3D0.0*ones(n,1); Y0(1)=3D0.0001;
   MaxTime=3D2910;
   m=3D0.001*ones(n,n); m=3Dm-diag(diag(m));
end

% Check all parameters are the right size
beta=3DCheckSize(beta,n,1,'beta');
gamma=3DCheckSize(gamma,n,1,'gamma');
nu=3DCheckSize(nu,n,1,'nu');
mu=3DCheckSize(mu,n,1,'mu');
m=3DCheckSize(m,n,n,'m');
X0=3DCheckSize(X0,n,1,'X0');
Y0=3DCheckSize(Y0,n,1,'Y0');

% Checks all the parameters are valid
CheckGreaterOrEqual(beta,0,'beta');
CheckGreaterOrEqual(gamma,0,'gamma');
CheckGreaterOrEqual(nu,0,'nu');
CheckGreaterOrEqual(mu,0,'mu');
CheckGreaterOrEqual(m,0,'m');
CheckGreaterOrEqual(X0,0,'X0');
CheckGreaterOrEqual(Y0,0,'Y0');

if sum(abs(diag(m)))&gt;0
    warning('Diagonal terms in movement matrix m are non-zero');
end

X=3DX0; Y=3DY0; All=3Dreshape([X';Y'],2*n,1);

% The main iteration=20
options =3D odeset('RelTol', 1e-4);
[t, pop]=3Dode45(@Diff_7_1,[0 MaxTime],All,options,n, beta, gamma, nu, mu,m=
);

X=3Dpop(:,1:2:(2*n)); Y=3Dpop(:,2:2:(2*n));

% plots the graphs with scaled colours
subplot(2,1,1)
h=3Dplot(t,X,'-g');
for i=3D1:n
    set(h(i),'Color',[0 0.5+0.5*(i-1)/(n-1) 0]);
end
xlabel 'Time';
ylabel 'Susceptibles'

subplot(2,1,2)=20
h=3Dplot(t,Y,'-r');
for i=3D1:n
    set(h(i),'Color',[0.5+0.5*(i-1)/(n-1) 0 0]);
end
legend(h);
xlabel 'Time';
ylabel 'Infectious'



% Calculates the differential rates used in the integration.
function dPop=3DDiff_7_1(t,pop, n, beta, gamma, nu, mu, m)

X=3Dpop(1:2:(2*n)); Y=3Dpop(2:2:(2*n));

dPop=3Dzeros(2*n,1);

% Note the different use of .* and *

dPop(1:2:(2*n))=3D nu -beta.*X.*Y - mu.*X + m*X - sum(m)'.*X;
dPop(2:2:(2*n))=3D beta.*X.*Y - gamma.*Y - mu.*Y + m*Y - sum(m)'.*Y;


% Does a simple check on the size and possible transpose
function Parameter=3DCheckSize(Parameter, L, W, str)

[l w]=3Dsize(Parameter);

if(l=3D=3DW &amp; w=3D=3DL)
    Parameter=3DParameter'; [l w]=3Dsize(Parameter);
end

if(l=3D=3D1 &amp; w=3D=3D1)
    warning('Parameter %s is a scaler value, expanding to size %d x %d',str=
,L,W);
    Parameter=3DParameter*ones(L,W);
    [l w]=3Dsize(Parameter);
end

if(l~=3DL | w~=3DW)
    error('Parameter %s is of size %d x %d and not %d x %d',str,l,w,L,W);
end


% Does a simple check on the value
function []=3DCheckGreaterOrEqual(Parameter, Value, str)

m=3Dfind(Parameter&lt;Value);
if length(m)&gt;0
    error('Parameter %s(%g) (=3D%g) is less than %g',str,m(1),Parameter(m(1=
)),Value);
end


</pre></body></html>
------MultipartBoundary--O7eOMlZSlJC40B9HEldG57t80vTSAZsSRTgVi59Jih------
