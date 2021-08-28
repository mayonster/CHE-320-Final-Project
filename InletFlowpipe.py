# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 01:37:39 2021

@author: ianna
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spop

R = 0.1225 #m
L = 0.0081 #km
dP = 18.43 #MJ/m^3
Vis0 = 16.355 #kJ*s/m^3
Gamma = -0.24
Lambda = 0.16
dPdz = dP/L #kJ/m^4

file1 = open("InletData.txt","a")#append mode

def Viscosity(dvdr):
    Factor1 = Lambda*dvdr
    Factor2 = Factor1*Factor1
    Vis = (1.+Factor2)**Gamma
    return Vis0*Vis

def Functions(v_Array,r_Array):
    f_Array = np.zeros_like(r_Array)
    f_Array[0] = v_Array[1]-v_Array[0]
    f_Array[-1] = v_Array[-1]
    for m in range(1,Bins):
        r_m = r_Array[m]
        v_mns = v_Array[m-1]
        v_pls = v_Array[m+1]
        v_m = v_Array[m]
        dvdr_m = 0.5*(v_pls-v_mns)/dr
        Vis_m = Viscosity(dvdr_m)
        dr_m = dr/r_m
        Pls = v_pls*(2.+dr_m)
        Mns = v_mns*(2.-dr_m)
        Extra = dPdzdrdr/Vis_m
        f_Array[m] = Pls+Mns-4.*v_m+Extra
    return f_Array


dr = 0.000025
dPdzdrdr = 2.*dPdz*dr*dr
r_Array = np.arange(0,R+dr,dr)
v_Array = np.zeros_like(r_Array)
v_Array = (dPdz/(4.*Vis0))*(R*R-r_Array*r_Array)
v_Array0 = (dPdz/(4.*Vis0))*(R*R-r_Array*r_Array)
Nodes = len(r_Array)
Bins = Nodes-1

Fs = Functions(v_Array,r_Array)
#print(Fs)


v_Array = spop.fsolve(Functions,v_Array,r_Array)

Fs = Functions(v_Array,r_Array)
#print(Fs)


dvdr_Array = np.gradient(v_Array,r_Array)
Viscosity_Array = Viscosity(dvdr_Array)
Tau_Array = dvdr_Array*Viscosity_Array
Force = Tau_Array[-1]
Force = Force*2*np.pi*R*L
print("Force",Force)
file1.write("Force: {}\n".format(Force))

Flowrate = np.trapz(v_Array*r_Array,r_Array)
Flowrate = 2.*np.pi*Flowrate
print("Flowrate: {}".format(Flowrate))
file1.write("Flowrate: {}\n".format(Flowrate))
Power = Flowrate*dP
print("Power",Power)
file1.write("Power: {}\n".format(Power))


plt.plot(r_Array,v_Array)
plt.axis([0,R,0,np.amax(v_Array)])
plt.xlabel("r [m]")
plt.ylabel("v [m/s]")
plt.show()


plt.plot(r_Array,Tau_Array)
plt.axis([0,R,np.amin(Tau_Array),0])
plt.xlabel("r [m]")
plt.ylabel("Tau [N/m\u00b2]")
plt.show()

print()
Line = "{:12s}\t{:12s}\t{:12s}\t\n".format("r [m]","vel [m/s]","Tau [N/m\u00b2]")
print(Line)
file1.write(Line)
for m in range(Nodes):
    Line = "{:12.6f}\t{:12.6f}\t{:12.4f}\t\n".format(r_Array[m],v_Array[m],Tau_Array[m])
    print(Line)
    file1.write(Line)


file1.close()