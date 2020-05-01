# Radiator Fins Toy Model

import numpy as np
import matplotlib.pyplot as plt
from pylab import *

def init_plotting():
	plt.rcParams['figure.figsize'] = (10,10)
	plt.rcParams['font.size'] = 10
	plt.rcParams['font.family'] = 'Times New Roman'
	plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
	plt.rcParams['axes.titlesize'] = 1.2*plt.rcParams['font.size']
	plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
	plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
	plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
	# plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
	plt.rcParams['xtick.major.size'] = 3    
	plt.rcParams['xtick.minor.size'] = 3
	plt.rcParams['xtick.major.width'] = 1
	plt.rcParams['xtick.minor.width'] = 1   
	plt.rcParams['ytick.major.size'] = 3
	plt.rcParams['ytick.minor.size'] = 3
	plt.rcParams['ytick.major.width'] = 1
	plt.rcParams['ytick.minor.width'] = 1
	plt.rcParams['legend.frameon'] = True
	plt.rcParams['legend.loc'] = 'best'
	plt.rcParams['axes.linewidth'] = 1

	plt.rcParams['lines.linewidth'] = 2.0 
	plt.rcParams['lines.markersize'] = 12

	plt.rcParams['figure.facecolor'] = 'white'
	plt.rcParams['axes.facecolor'] = 'white'
	#plt.rcParams['axes.color_cycle'] = ['b', 'r', 'g','pink','orange','darkgreen','purple']

	plt.rcParams['grid.color'] = 'k'
	plt.rcParams['grid.linestyle'] = ':'
	plt.rcParams['grid.linewidth'] = 0.5

	#plt.gca().spines['right'].set_color('None')
	#plt.gca().spines['top'].set_color('None')
	plt.gca().xaxis.set_ticks_position('bottom')
	plt.gca().yaxis.set_ticks_position('left')
init_plotting()

# Fah=0.
# Foh=0.
# Cs=0.
# Cl=0.

# nlayers=2
# psurf=1000.
# p=np.linspace(psurf,0.,nlayers)
# h=7.
# gamma=6.2
# zxT=16.
# sb=5.67e-8

# Fsens=0.
# Idown[0]=393.
# r=0.8
# alpha=0.8
# Sc=370
# ustar=4
# Cd=0.002

# # (2.1)
# H=Fv[nlayers]+Fah+Foh 	# H		=	rate of energy accumulation by column
# 					# Fv0	=	net vertical flux out of BOA
# 					# Fah	=	net horizontal atmospheric heat transport (latent+sensible) into column (per unit cross-section area)
# 					# Foh	=	net horizontal oceanic heat transport into column (per unit cross-section area)

# # (2.2a)
# S=Sc+Cs 	# S 	=	absorbed solar radiation measured at TOA
# 			# Sc 	= 	estimate of absorbed solar radiation without clouds
# 			# Cs 	=	cloud shortwave forcing
# # (2.2b)
# Iup[nlayers]=Iupc[nlayers]-Cl 	# Iup[nlayers] 	=	OLR
# 								# Iupc[nlayers]	= 	estimate of OLR without clouds
# 								# Cl 			=	cloud longwave forcing

# # (2.3)
# Fvinf=Sc-Iupc[nlayers]+(Cs-Cl) 	# Fvinf 	= 	net vertical flux into TOA 

# # (2.4)
# zx=-h*log(p/psurf)
# for i in range(nlayers):
# 	if(zx[i]<zxT):
# 		T[i]=T[0]-gamma*zx
# 	else:
# 		T[i]=T[0]-gamma*zxT

# Ts=T[0]+1.

# # (2.5a)
# Fv[0]=-Foh	# Fv	=	net vertical flux

# # (2.5b)
# Fv[0]=alpha*Sc+Cstars+(Idown[0]-sb*Ts**4)-E-Fsens 	# alpha 	=	SW absorption coefficient representing portion of clear sky TOA solar flux that reaches the ground (different to albedo)
# 														# Cstars	=	effect of clouds on surface insolation (could be different to Cs)
# 														# Idown 	=	downwelling IR flux
# 														# E 		=	evaporative heat flux
# 														# Fsens 	=	sensible heat flux	

# # (2.6a)
# E=rho*L*Cd*ustar*qstar	# rho	= 	air density in boundary layer
# 						# L 	=	latemt heat of vaporisation
# 						# Cd 	=	drag coefficient
# 						# ustar = 	characteristic velocity fluctuation
# 						# qstar = 	characteristic scale of fluctuation of water vapour mass mixing ratio

# # (2.6b)
# qstar=qsat(Ts)-r*qsat(T[0])

# # (2.7)
# if(Ts>300):
# 	Cstars=-a*(Ts-300)
# else:
# 	Cstars=0.


# SEB=alpha*Sc+Cstars+(Idown[0]-sb*Ts**4.)

# # (3.1)
# Qv=Fv[nlayers]-Fv[0]	# Qv 	=	energy added to atmospheric column due to convergence of vertical flux

# # (3.2)
# H=Qv+Fah+(Fv[0]+Foh)

# # (3.3)
# H=Gv+Fah

# # (3.4)
# Qv=Sc-Iupc[nlayers]+(Cs-Cl)+Foh

# # (3.5)
# Fah=Faexp+meanQv-Qv

# # (3.6)
# H=meanQv+Faexp

# # (3.7a)
# Q=Qv+(P-E) 	# P 	=	latent heat release due to precip
# 			# E 	=	evaporative heat flux

# # (3.7b)
# Fah=(P-E)+Fahs 	# Fahs 	=	atmospheric sensible heat transport

########################################################################################################################################################################################################

def qsat(T,p):
	Tc=T-273.15
	es=6.1094*np.exp(17.625*Tc/(Tc+243.04))
	qsat=es/9
	return qsat

nlayers=10
psurf=1000.
h=7.
gamma=6.2
zxT=16.
Faexp=0.
Foh1=0.
Foh2=0.
A1=1.
A2=2.
Sc=370.
Cs=0.
Cl=0.
Cp=1006. # eng toolbox
r2=0.7
L=2264.705 # wiki
deltheta=85.
alpha=0.8
Cd=0.002
e2=0.25
Cstars=0
Fsens=0.
Foh=0

zxrad2=h/2.
sb=5.67e-8

# Fah1=(A2/A1)*Qv2

T1=np.zeros(nlayers)
T1[0]=280.
psurf=1000.
p=np.linspace(psurf,1.,nlayers)
zx=-h*np.log(p/psurf)
for i in range(1,nlayers):
	if(zx[i]<zxT):
		T1[i]=T1[0]-gamma*zx[i]
	else:
		T1[i]=T1[0]-gamma*zxT

T2=T1
r1=0.75
data=np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/x[T0]_y[OLR]_r75_Pierrehumbert95.txt',delimiter=',')
xp=data[:,0]
yp=data[:,1]
OLR1=np.interp(T1[0],xp,yp)
Ts1=T1[0]+1.
Qv1=Sc-OLR1+(Cs-Cl)+Foh1
Tbar=T2[np.argmin(abs(zx-zxrad2))]
Ts2=Ts1
for i in range(30):
	if(i==0):
		qstar2=r2*qsat(Ts2,psurf)
	b=L*qstar2/(Cp*deltheta)
	Q2=(1-alpha)*Sc+e2*sb*Ts2**4-2.*e2*sb*Tbar**4
	E2=-b*(Q2+Faexp)
	Idown0=393. # from p1788, not proper
	# Idown0=e2*sb*Tbar**4
	Fv0=alpha*Sc+Cstars+(Idown0-sb*Ts2**4)-E2-Fsens 
	SEB2=Fv0+Foh
	Ts2=((1./(1.-e2*b)*((1.+(e2-1)*alpha)*Sc+e2*(e2-2.)*sb*Tbar**4+Faexp+e2*Foh2))/sb)**0.25
	Qv2=( 1. - b ) / ( 1. - e2 * b ) * ( ( 1. + ( e2 - 1 ) * alpha ) * Sc + e2 * ( e2 - 2. ) * sb * Tbar**4 + Faexp + e2 * Foh2 )
	Fah1=(A2/A1)*Qv2	
	OLR1=np.interp(T1[0],xp,yp)
	plt.figure(1)
	plt.plot(T1[0],OLR1,'o',c='b')
	plt.plot(T1[0],Sc+Fah1,'o',c='r')
	T1[0]+=(Sc+Fah1-OLR1)*0.1
	Ts1=T1[0]+1
	for i in range(1,nlayers):
		if(zx[i]<zxT):
			T1[i]=T1[0]-gamma*zx[i]
		else:
			T1[i]=T1[0]-gamma*zxT
	T2=T1
	Tbar=T2[np.argmin(abs(zx-zxrad2))]

################################################################################################################################################################
show()