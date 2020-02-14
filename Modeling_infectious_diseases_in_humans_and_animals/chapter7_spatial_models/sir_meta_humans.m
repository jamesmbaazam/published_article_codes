From: <Saved by Blink>
Snapshot-Content-Location: http://www.modelinginfectiousdiseases.org/
Subject: Modeling Infectious Diseases in Humans and Animals, by Keeling & Rohani
Date: Wed, 12 Feb 2020 14:34:37 -0000
MIME-Version: 1.0
Content-Type: multipart/related;
	type="text/html";
	boundary="----MultipartBoundary--jdk7jmVQ2OCSADMIC4x0VY6gUr0vsCejm0Z5Q8dQsX----"


------MultipartBoundary--jdk7jmVQ2OCSADMIC4x0VY6gUr0vsCejm0Z5Q8dQsX----
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
------MultipartBoundary--jdk7jmVQ2OCSADMIC4x0VY6gUr0vsCejm0Z5Q8dQsX----
Content-Type: text/html
Content-ID: <frame-A68C0874F133FB69827CA7B11B8F8D46@mhtml.blink>
Content-Transfer-Encoding: quoted-printable
Content-Location: http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter7/Program_7.2/Program_7_2.m

<html><head><meta http-equiv=3D"Content-Type" content=3D"text/html; charset=
=3Dwindows-1252"></head><body><pre style=3D"word-wrap: break-word; white-sp=
ace: pre-wrap;">function [t,X,Y] =3D Program_7_2(n,beta,gamma,l,r,N0,X0,Y0,=
MaxTime)
%
%=20
%
% Program_7_2( n,beta, gamma, l, r, N0, X0, Y0, MaxTime)
%      This is the MATLAB version of program 7.2 from page 242 of=20
% "Modeling Infectious Disease in humans and animals"=20
% by Keeling &amp; Rohani.
%=20
% It is the SIR epidemic in a metapopulationFor simplicity births and death=
s have been=20
% ignored, and we work with numbers of individuals.
% Y[i][j] refers to infected individual who are currently in i but live in =
j..
%

% Sets up default parameters if necessary.
if nargin =3D=3D 0
   n=3D5;
   beta=3D1.0*ones(n,1);
   gamma=3D0.3*ones(n,1);
   N0=3D1000*ones(n,1);
   X0=3D800*ones(n,1);
   Y0=3D0.0*ones(n,1); Y0(1)=3D1;
   MaxTime=3D60;
   l=3D0.1*(diag(ones(1,n-1),1)+diag(ones(1,n-1),-1));
   r=3D2*ones(n,n); r=3Dr-diag(diag(r));
end

% Check all parameters are the right size
beta=3DCheckSize(beta,n,1,'beta');
gamma=3DCheckSize(gamma,n,1,'gamma');
l=3DCheckSize(l,n,n,'m');
r=3DCheckSize(r,n,n,'m');
N0=3DCheckSize(N0,n,1,'S0');
X0=3DCheckSize(X0,n,1,'S0');
Y0=3DCheckSize(Y0,n,1,'I0');

% Checks all the parameters are valid
CheckGreaterOrEqual(beta,0,'beta');
CheckGreaterOrEqual(gamma,0,'gamma');
CheckGreaterOrEqual(l,0,'m');
CheckGreaterOrEqual(r,0,'m');
CheckGreaterOrEqual(N0,0,'S0');
CheckGreaterOrEqual(X0,0,'S0');
CheckGreaterOrEqual(Y0,0,'I0');

if sum(abs(diag(l)))&gt;0
    warning('Diagonal terms in movement matrix l are non-zero');
end
if sum(abs(diag(r)))&gt;0
    warning('Diagonal terms in return rate matrix r are non-zero');
end
m=3Dfind(X0+Y0&gt;N0);
if length(m)&gt;0
    error('Population %g, susceptibles + infecteds are greater than populat=
ion size',m(1));
end

X=3Dzeros(n,n); Y=3Dzeros(n,n); NN=3Dzeros(n,n);
X=3Ddiag(X0); Y=3Ddiag(Y0); NN=3Ddiag(N0); All=3Dreshape([X Y NN],3*n*n,1);

% The main iteration=20
options =3D odeset('RelTol', 1e-4);
[t, pop]=3Dode45(@Diff_7_1,[0 MaxTime],All,options,n, beta, gamma, l, r);

X=3Dzeros(length(t),n); Y=3DX;
for i=3D1:n
    X(:,i)=3Dsum(pop(:,i:n:(n*n))')';
    Y(:,i)=3Dsum(pop(:,(n*n)+[i:n:(n*n)])')';
end

% plots the graphs with scaled colours
subplot(2,1,1)
h=3Dplot(t,X,'-g');
for i=3D1:n
    set(h(i),'Color',[0 0.5+0.5*(i-1)/(n-1) 0]);
end
legend(h);
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
function dPop=3DDiff_7_1(t,pop, n, beta, gamma, l, r)

X=3Dreshape(pop(1:(n*n)),n,n); Y=3Dreshape(pop((n*n)+[1:(n*n)]),n,n);=20
NN=3Dreshape(pop(2*(n*n)+[1:(n*n)]),n,n);

dX=3Dzeros(n,n); dY=3Dzeros(n,n); dNN=3Dzeros(n,n);

% Note the different use of .* and *

sumY=3Dsum(Y')'; sumNN=3Dsum(NN')';

% First do the off diagonals.
dX =3D -X.*((beta.*sumY./sumNN)*ones(1,n)) + l.*(ones(n,1)*diag(X)')   - r.=
*X;
dY =3D X.*((beta.*sumY./sumNN)*ones(1,n)) - Y.*(gamma*ones(1,n)) + l.*(ones=
(n,1)*diag(Y)')   - r.*Y;
dNN =3D  l.*(ones(n,1)*diag(NN)')   - r.*NN;

% Now do diagonals
DX =3D -diag(X).*(beta.*sumY./sumNN) - sum(l)'.*diag(X) + sum(r.*X)';
DY =3D diag(X).*(beta.*sumY./sumNN) - diag(Y).*gamma - sum(l)'.*diag(Y) + s=
um(r.*Y)';
DNN =3D - sum(l)'.*diag(NN) + sum(r.*NN)';

dX =3D dX - diag(diag(dX)) + diag(DX);
dY =3D dY - diag(diag(dY)) + diag(DY);
dNN =3D dNN - diag(diag(dNN)) + diag(DNN);

dPop=3Dreshape([dX dY dNN],3*n*n,1);

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
------MultipartBoundary--jdk7jmVQ2OCSADMIC4x0VY6gUr0vsCejm0Z5Q8dQsX------
