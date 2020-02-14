From: <Saved by Blink>
Snapshot-Content-Location: http://www.modelinginfectiousdiseases.org/
Subject: Modeling Infectious Diseases in Humans and Animals, by Keeling & Rohani
Date: Wed, 12 Feb 2020 14:34:22 -0000
MIME-Version: 1.0
Content-Type: multipart/related;
	type="text/html";
	boundary="----MultipartBoundary--fzUl63nJXXRS9NakHpmVrVRPNEvzn9lN2G1NfZW8vX----"


------MultipartBoundary--fzUl63nJXXRS9NakHpmVrVRPNEvzn9lN2G1NfZW8vX----
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
------MultipartBoundary--fzUl63nJXXRS9NakHpmVrVRPNEvzn9lN2G1NfZW8vX----
Content-Type: text/html
Content-ID: <frame-A68C0874F133FB69827CA7B11B8F8D46@mhtml.blink>
Content-Transfer-Encoding: quoted-printable
Content-Location: http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter7/Program_7.2/Program_7_2.py

<html><head><meta http-equiv=3D"Content-Type" content=3D"text/html; charset=
=3Dwindows-1252"></head><body><pre style=3D"word-wrap: break-word; white-sp=
ace: pre-wrap;">#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 7.2 from page 242 of  #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling &amp; Rohani.										   #
###																   #
### It is the SIR epidemic in a metapopulationFor simplicity births#
### and deaths have been ignored, and we work with numbers of      #
### individuals.                                                   #
### Y[i][j] refers to infected individual who are currently in i   #
### but live in j..                                                #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl
from matplotlib.font_manager import FontProperties as fmp

n=3D5
beta=3D1.0*np.ones(n);
gamma=3D0.3*np.ones(n);
N0=3Dnp.zeros(n*n);
X0=3Dnp.zeros(n*n);
for i in np.arange(0,n*n,n+1):
	N0[i]=3D1000.0
	X0[i]=3D800.0

Y0=3Dnp.zeros(n*n); Y0[0]=3D1.0;
ND=3DMaxTime=3D60.
TS=3D1.0

l=3Dnp.zeros((n,n));r=3Dnp.zeros((n,n))
for i in range(n):
	for j in range(n):
		if abs(i-j)=3D=3D1:=20
			l[i][j]=3D0.1
r=3D2*np.ones((n,n)); r=3Dr-np.diag(np.diag(r));

INPUT0=3Dnp.hstack((X0,Y0,N0))
INPUT =3D np.zeros((3*n*n))
for i in range(n*n):
	INPUT[3*i]=3DINPUT0[i]
	INPUT[1+3*i]=3DINPUT0[n*n+i]
	INPUT[2+3*i]=3DINPUT0[2*n*n+i]

def diff_eqs(INP,t): =20
	'''The main set of equations'''
	Y=3Dnp.zeros((3*n*n))
	V =3D INP  =20
	sumY=3Dnp.zeros(n)
	sumN=3Dnp.zeros(n)
=09
	## Calculate number currently in Subpopulation i
	for i in range(n):
		sumY[i]=3D0.0;sumN[i]=3D0.0
		for j in range(n):
			k=3D3*(j+i*n);
			sumN[i]+=3DV[2+k];
			sumY[i]+=3DV[1+k];=09
		=09
	## Set all rates to zeros
	for i in range(n):
		for j in range(n):
			k=3D3*(j+i*n);
			Y[k]=3D0; Y[1+k]=3D0; Y[2+k]=3D0
=09
	for i in range(n):
		for j in range(n):	=09
			## Calculate the rates
			k =3D 3 * (j+i*n)=20
			K =3D 3 * (i+j*n)
			h =3D 3 * (i+i*n)
			H =3D 3 * (j+j*n)
		=09
			Y[k] -=3D (beta[i]*V[k]*(sumY[i]/sumN[i]))
			Y[k+1] +=3D (beta[i]*V[k]*(sumY[i]/sumN[i]))
			Y[k+1] -=3D (gamma[i]*V[k+1])
		=09
			## Movement
			Y[h] +=3D r[j][i]*V[K]
			Y[h] -=3D l[j][i]*V[h]
		=09
			Y[h+1] +=3D r[j][i]*V[K+1]
			Y[h+1] -=3D l[j][i]*V[h+1]
		=09
			Y[h+2] +=3D r[j][i]*V[K+2]
			Y[h+2] -=3D l[j][i]*V[h+2]
		=09
			Y[k] +=3D l[i][j]*V[H]
			Y[k] -=3D r[i][j]*V[k]
		=09
			Y[1+k] +=3D l[i][j]*V[1+H]
			Y[1+k] -=3D r[i][j]*V[1+k]
		=09
			Y[2+k] +=3D l[i][j]*V[2+H]
			Y[2+k] -=3D r[i][j]*V[2+k]
	return Y   # For odeint

t_start =3D 0.0; t_end =3D ND; t_inc =3D TS
t_range =3D np.arange(t_start, t_end+t_inc, t_inc)
t_course =3D spi.odeint(diff_eqs,INPUT,t_range)
tc =3D t_course

### Plotting
totalS=3Dnp.zeros((len(tc),5))
totalI=3Dnp.zeros((len(tc),5))

for i in range(n):
	for j in range(n):
		k=3D3*(j+i*n);
		totalS[:,i]+=3Dtc[:,k]
		totalI[:,i]+=3Dtc[:,k+1]


#print len(totalS)
pl.subplot(211)
for i in range(5):
	pl.plot(t_range,totalS[:,i], label=3D('data %s' %(i+1)), color=3D(0.3,i/10=
.+0.5,0.1))
pl.xlabel('Time')
pl.ylabel('Susceptibles')
pl.legend(loc=3D1,prop =3D fmp(size=3D'smaller'))
pl.subplot(212)
for i in range(5):
	pl.plot(t_range,totalI[:,i], label=3D('data %s' %(i+1)), color=3D(0.8,i/10=
.+0.,0.3))
pl.xlabel('Time')
pl.ylabel('Infectious')
pl.legend(loc=3D1,prop =3D fmp(size=3D'smaller'))

pl.show()
</pre></body></html>
------MultipartBoundary--fzUl63nJXXRS9NakHpmVrVRPNEvzn9lN2G1NfZW8vX------
