#TEO_RV_Model.py
#This code overplots RV data with an RV model.
#The calculations follow Perryman 2011, The Exoplanet Handbook, Cambridge University Press, pp. 9-13.

import math
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#Import RV data to be fit. Should be of the form []
Data = np.loadtxt("RVDataGood.txt")

#Constants (MKS units)
G = 6.67408e-11
MJup = 1.89813e27 
MSun = 1.98855e30 #https://sites.google.com/site/mamajeksstarnotes/basic-astronomical-data-for-the-sun

#Independently Measured System Parameters
Ms = 1.364 #Host star's mass (M_Sun)

#Fit parameters
n = 30.0

Pmin = 3700.0      #orbital period (d)
Pmax = 4000.0
dP = (Pmax-Pmin)/n

tpmin = 55400.0    #time of periastron passage (d) (JD-24000000)
tpmax = 55700.0
dtp = (tpmax-tpmin)/n

emin = 0.2   #eccentricity (unitless)
emax = 0.4
de = (emax-emin)/n

Mpsinimin = 0.0   #planet's mass attenuated by sini (M_Jup)
Mpsinimax = 1.0
dMpsini = (Mpsinimax-Mpsinimin)/n

omegamin = 200.0   #longitude of periastron (degrees) (angle from the ascending node to periastron)
omegamax = 300.0
domega = (omegamax-omegamin)/n

gammamin = 0.0    #systematic velocity (m/s) (baseline offset)
gammamax = 2.0
dgamma = (gammamax-gammamin)/n

dmax = 0.0         #linear trend parameter (m/(s*d)) (baseline slope)
dmin = 0.0 
dd = (dmax-dmin)/n

#Initialize the fit parameters
P = Pmin
tp = tpmin
e = emin
Mpsini = Mpsinimin
omega = omegamin
gamma = gammamin
d = dmin
chi2nu = 0.0
Fits = np.array([[Pmin,tpmin,emin,Mpsinimin,omegamin,gammamin,dmin,chi2nu]])

#Compute model and chi2
P = Pmin
while P<=Pmax:
    tp = tpmin
    while tp<=tpmax:
        e = emin
        while e<=emax:
            Mpsini = Mpsinimin
            while Mpsini<=Mpsinimax:
                #while omega<=omegamax:
                    #while gamma<=gammamax:
                        #while d<=dmax:
                i = 0
                chi2 = 0.0
                for t in Data[:,0]:
                    assini = (Mpsini*MJup)*(G**0.5*(P*24.0*3600.0)/(2*np.pi*(Ms*MSun)))**(2.0/3.0)  #star's orbital semi-major axis projected into the plane of the sky (MKS) 
                    K = 2*np.pi*assini/((P*24*3600)*(1-e**2)**(0.5))  #K radial velocity semi-amplitude (km/s)
                    MA = (2*np.pi/P)*(t-tp)  #mean anomoly 
                    def fnc(x): 
                        return x-e*np.sin(x)-MA
                    E = fsolve(fnc,(np.pi/2))  #eccentric anomoly
                    v=2.0*math.atan2(math.sqrt(1+e)*math.tan(E/2),math.sqrt(1-e))  #true anomoly (angle between periastron and the orbiting body at a given time)
                    vr = K*(np.cos(omega*np.pi/180+v)+e*np.cos(omega*np.pi/180))+gamma+d*(t-tp)  #radial velocity (m/s)
                    chi2 = chi2+((vr-Data[i,1])/Data[i,2])**2
                    i = i+1
                    chi2nu = chi2/(len(Data)-7)
                Fits = np.append(Fits,[[P,tp,e,Mpsini,omega,gamma,d,chi2nu]],axis=0)
                            #d = d+dd
                        #gamma = gamma+dgamma
                    #omega = omega+domega
                Mpsini = Mpsini+dMpsini
            e = e+de
        tp = tp+dtp
    P = P+dP

Fits = Fits[1:,:]
print Fits
Row = np.argmin(Fits[:,7])
print Row
BestFit = Fits[Row]
print BestFit

#Overplot the best fit model and data

#Best fit parameters
P = BestFit[0]
tp = BestFit[1]
e = BestFit[2]
Mpsini = BestFit[3]
omega = BestFit[4]
gamma = 2
d = 0

#Non-free parameters
assini = (Mpsini*MJup)*(G**0.5*(P*24.0*3600.0)/(2*np.pi*(Ms*MSun)))**(2.0/3.0)
K = 2*np.pi*assini/((P*24*3600)*(1-e**2)**(0.5))
ap = assini*(Ms*MSun)/(Mpsini*MJup)

#Generate the time points for the model plot
margin = 500      #desired extension of the model both before and after the data (d)
Ntimesteps = 500  #desired number of time points for the RV model plot
NData = len(Data)-1
tmin = Data[0,0]-margin
tmax = Data[NData,0]+margin
timestep = (tmax-tmin)/Ntimesteps
Time = np.arange(tmin,tmax,timestep)
RVBest = np.array([])

#Generate the RV points for the model plot
for t in Time:
    MA = (2*np.pi/P)*(t-tp)  #mean anomoly
    def fnc(x):
        return x-e*np.sin(x)-MA
    E = fsolve(fnc,(np.pi/2))  #eccentric anomoly
    #v = np.arccos((np.cos(E)-e)/(1-e*np.cos(E)))
    v=2.0*math.atan2(math.sqrt(1+e)*math.tan(E/2),math.sqrt(1-e))  #true anomoly (angle between periastron and the orbiting body at a given time)
    vr = K*(np.cos(omega*np.pi/180+v)+e*np.cos(omega*np.pi/180))+gamma+d*(t-tp)  #radial velocity (m/s)
    RVBest = np.append(RVBest,vr)

#Generate the plot
plt.plot(Time, RVBest, "r")
plt.plot(Data[:,0], Data[:,1], ".k")
plt.errorbar(Data[:,0], Data[:,1], yerr=Data[:,2], linestyle="none", color="k", capsize=2.5)
plt.show