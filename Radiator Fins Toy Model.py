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

def qsat(T,p):
	Tc=T-273.15
	es=6.1094*np.exp(17.625*Tc/(Tc+243.04))
	qsat=621.97*(es)/(p-es)*1e-3
	return qsat

nlayers=100
psurf=1000.
h=7.
gamma=6.2
zxT=16.
Faexp=0.
Foh1=0.
Foh2=0.
A1=1.
A2=1.
Sc=370.
Cs=0.
S=Sc+Cs
Cl=0.
Cp=1.006 # eng toolbox
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
sb=5.670374419e-8

data=np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/x[T0]_y[OLR]_r75_Pierrehumbert95.txt',delimiter=',')
# data=np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/x[T0]_y[OLR]_r100_Pierrehumbert95.txt',delimiter=',')
xp=data[:,0]
yp=data[:,1]

# Fah1=(A2/A1)*Qv2

i_leg=0

# e2s=np.linspace(0.35,0.45,100)
e2s=np.linspace(0.1,1.,10)
# print e2s
# e2s=[0.1]
for e2 in e2s:
	T=np.zeros(nlayers)
	T[0]=300.
	psurf=1000.
	p=np.linspace(psurf,1.,nlayers)
	zx=-h*np.log(p/psurf)
	for i in range(1,nlayers):
		if(zx[i]<zxT):
			T[i]=T[0]-gamma*zx[i]
		else:
			T[i]=T[0]-gamma*zxT
		Tbar=T[np.argmin(abs(zx-zxrad2))]
	Ts1=T[0]+1.
	Ts2=T[0]+1.
	qstar2=0
	b=0
	E2=0
	SEB2=0
	Qv2=0
	Fah1=0
	Fah2=0
	Qv1=0
	H1=0
	H2=0
	OLR1=0
	OLR2=0
	TOA_budg_1=0
	TOA_budg_2=0
	L=2500.
	# Qv1=Sc-OLR1+(Cs-Cl)+Foh1
	# Tbar=T2[np.argmin(abs(zx-zxrad2))]
	# Ts2=Ts1
	# qstar2=r2*qsat(Ts2,psurf)
	# T=np.zeros(nlayers)
	# t0s=np.linspace(250,350,20)
	# T[0]=300.

	# for t0 in t0s:
		# T[0]=t0
	for i_time in range(10000):
		# if(i<5):
		Tc=T[0]-273.15
		L=2500.8-2.36*Tc+0.0016*Tc**2-0.00006*Tc**3.
		qstar2=r2*qsat(Ts2,psurf)
		b=(L*qstar2/(Cp*deltheta))
		Q2=(1.-alpha)*Sc+e2*sb*Ts2**4.-2.*e2*sb*Tbar**4.
		E2=-b*(Q2+Faexp)
		# Idown0=393. # from p1788, not proper
		# Idown0=e2*sb*Tbar**4
		# Fv0=alpha*Sc+Cstars+(Idown0-sb*Ts2**4)-E2-Fsens 
		SEB2=(1.)/(1.-e2*b) * ((alpha+(1.-alpha)*b)*Sc+(1.-2.*b)*e2*sb*Tbar**4.+b*Faexp+Foh2)
		# Ts2=(((1.)/(1.-e2*b) * ((alpha+(1.-alpha)*b)*Sc+(1.-2.*b)*e2*sb*Tbar**4.+b*Faexp+Foh2))/sb)**0.25
		# dTs2=np.clip(SEB2,-0.01,0.01)
		Qv2=( 1. - b ) / ( 1. - e2 * b ) * ( ( 1. + ( e2 - 1. ) * alpha ) * Sc + e2 * ( e2 - 2. ) * sb * Tbar**4. + Faexp + e2 * Foh2 ) - Faexp
		Fah1=(A2/A1)*Qv2	
		Qv1=(-(Qv2+Faexp)*(A2)/(A1+A2))*(A2)/(A1+A2)-Faexp
		Fah2=-Fah1
		H1=(Qv1+Faexp)*(A1/(A1+A2))
		H2=(Qv2+Faexp)*(A2/(A1+A2))
		# T[0]+=H2*0.1
		OLR1=np.interp(T[0],xp,yp)
		OLR2=e2*sb*Tbar**4.
		# TOA_budg_1=Sc+Fah1-OLR1
		TOA_budg_1=Sc+Fah1-OLR1
		TOA_budg_2=S+Fah2-OLR2
		dT0=np.clip(TOA_budg_1*0.001,-1e3,1e3)
		T[0]+=dT0
		# T[0]+=TOA_budg_1*0.01
		# T[0]+=(TOA_budg_1+TOA_budg_2)*0.01
		for i in range(1,nlayers):
			if(zx[i]<zxT):
				T[i]=T[0]-gamma*zx[i]
			else:
				T[i]=T[0]-gamma*zxT
		Ts1=T[0]+1.
		# Ts2+=SEB2*1e-3
		# Ts2=(((1.)/(1.-e2*b) * ((alpha+(1.-alpha)*b)*Sc+(1.-2.*b)*e2*sb*Tbar**4.+b*Faexp+Foh2))/sb)**0.25
		Ts2+=((((1.)/(1.-e2*b) * ((alpha+(1.-alpha)*b)*Sc+(1.-2.*b)*e2*sb*Tbar**4.+b*Faexp+Foh2))/sb)**0.25 - Ts2)*1e-3

		# Ts2-=E2*0.01
		Tbar=T[np.argmin(abs(zx-zxrad2))]
		Tc=T[0]-273.15
		L=2500.8-2.36*Tc+0.0016*Tc**2-0.00006*Tc**3.
		if(i_time%50==0):
			print '{:4.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f}'.format(b,Tbar,Ts1,Ts2,TOA_budg_1,TOA_budg_2,SEB2)
		# plt.figure(1)
		# if(b<0.5):
		# 	plt.subplot(121)
		# 	plt.plot(Tbar,Ts2,'o',c='b')
		# else:
		# 	plt.subplot(122)
		# 	plt.plot(Tbar,Ts2,'o',c='b')
		# plt.plot(T[0],H1,'o',c='r')
		# plt.plot(T[0],H2,'o',c='b')
		# plt.axhline(0)
		# print H1, H2, T[0]
		# plt.figure(1)
		# plt.plot(T[0],OLR1,'o',c='b')
		# plt.plot(T[0],Sc+Fah1,'o',c='r')
		# if(abs(H1)<1e-6 and abs(H2)<1e-6):
		# if(abs(TOA_budg_1)<1e-1 and abs(SEB2)<1e-1):
		if(abs(TOA_budg_1)<1e-1):
			print 'eqb reached'
			plt.figure(1)
			plt.plot(e2,Ts1,'o',c='r',label='Ts1 model')
			plt.plot(e2,Ts2,'o',c='b',label='Ts2 model')
			plt.xlabel('e2')
			plt.ylabel('T')
			if(i_leg==0):
				plt.legend()
			i_leg+=1
			# plt.plot(T'bar,Ts2,'o')
			# plt.plot(T[0],OLR1,'o',c='b')
			# plt.plot(T[0],Sc+Fah1,'o',c='r')
			print b, Ts1, Ts2,E2 #put area weighted mean T here.
			break
	print '-----------------------------------------'

data=np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/x[e2]_y[Ts1].txt',delimiter=',')
xp=data[:,0]
yp=data[:,1]
plt.figure(1)
plt.plot(xp,yp,c='r')
data=np.genfromtxt('/Users/nickedkins/Dropbox/Figure Data/x[e2]_y[Ts2].txt',delimiter=',')
xp=data[:,0]
yp=data[:,1]
plt.figure(1)
plt.plot(xp,yp,c='b')

		
		# print H1, H2,TOA_budg_1, TOA_budg_2
		# T1[0]+=TOA_budg_1*0.001
		# T2[0]+=TOA_budg_2*0.001
		# Ts1=T1[0]+1
		# plt.figure(1)
		# plt.plot(Ts1,Ts2,'o')
		# for i in range(1,nlayers):
		# 	if(zx[i]<zxT):
		# 		T1[i]=T1[0]-gamma*zx[i]
		# 	else:
		# 		T1[i]=T1[0]-gamma*zxT
		# T2=T1
		# Tbar=T2[np.argmin(abs(zx-zxrad2))]
		# Tbar=T[np.argmin(abs(zx-zxrad2))]
		# if(abs(TOA_budg_1)<1e-3):
		# if(abs(H1)<1.0 and abs(H2)<1.0):
		# 	print 'eqb'
		# 	# plt.figure(1)
		# 	# # if(b>0.5):
		# 	# # 	plt.plot(Ts2,E2,'o')
		# 	# plt.plot(Ts2,b,'o')
		# 	plt.plot(e2,Ts1,'o',c='r')
		# 	plt.plot(e2,Ts2,'o',c='b')
		# 	break

################################################################################################################################################################
show()