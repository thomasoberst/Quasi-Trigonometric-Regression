#TEO_RV_Model.py
#This code overplots RV data with an RV model.
#The calculations follow Perryman 2011, The Exoplanet Handbook, Cambridge University Press, pp. 9-13.

import math
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#Import RV data
Data = np.loadtxt("RVDataGood.txt")

#Constants (MKS units)
G = 6.67408e-11
MJup = 1.89813e27 
MSun = 1.98855e30 #https://sites.google.com/site/mamajeksstarnotes/basic-astronomical-data-for-the-sun

#Independently Measured System Parameters
Ms = 1.364 #Host star's mass (M_Sun)

#Free parameters
P = 3700.0      #orbital period (d)
tp = 55700.0    #time of periastron passage (d) (JD-24000000)
e = 0.99        #eccentricity (unitless)
Mpsini = 0.1   #planet's mass attenuated by sini (M_Jup)
omega = 0.0   #longitude of periastron (degrees) (angle from the ascending node to periastron)
gamma = 0.0    #systematic velocity (m/s) (baseline offset)
d = 0.0         #linear trend parameter (m/(s*d)) (baseline slope)

#Non-free parameters
assini = (Mpsini*MJup)*(G**0.5*(P*24.0*3600.0)/(2*np.pi*(Ms*MSun)))**(2.0/3.0)  #star's orbital semi-major axis projected into the plane of the sky (MKS)
K = 2*np.pi*assini/((P*24*3600)*(1-e**2)**(0.5))                                #K radial velocity semi-amplitude (km/s)
ap = assini*(Ms*MSun)/(Mpsini*MJup)                                             #planet's orbital semimajor axis (m)

#Generate the time points for the model plot
margin = 500      #desired extension of the model both before and after the data (d)
Ntimesteps = 500  #desired number of time points for the RV model plot
NData = len(Data)-1
tmin = Data[0,0]-margin
tmax = Data[NData,0]+margin
timestep = (tmax-tmin)/Ntimesteps
Time = np.arange(tmin,tmax,timestep)
RV = np.array([])

#Generate the RV points for the model plot
for t in Time:
    MA = (2*np.pi/P)*(t-tp)  #mean anomoly
    def fnc(x):
        return x-e*np.sin(x)-MA
    E = fsolve(fnc,(np.pi/2))  #eccentric anomoly
    #v = np.arccos((np.cos(E)-e)/(1-e*np.cos(E)))
    v=2.0*math.atan2(math.sqrt(1+e)*math.tan(E/2),math.sqrt(1-e))  #true anomoly (angle between periastron and the orbiting body at a given time)
    vr = K*(np.cos(omega*np.pi/180+v)+e*np.cos(omega*np.pi/180))+gamma+d*(t-tp)  #radial velocity (m/s)
    RV = np.append(RV,vr)

#Generate the plot
plt.plot(Time, RV, "r")
plt.plot(Data[:,0], Data[:,1], ".k")
plt.errorbar(Data[:,0], Data[:,1], yerr=Data[:,2], linestyle="none", color="k", capsize=2.5)
plt.show