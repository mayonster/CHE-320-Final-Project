# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 18:34:47 2021

@author: ianna
"""

import numpy as np
import matplotlib.pyplot as plt
#import time

k = 0.0751 # [/s] 
alpha = 0.21
beta = 1.36

#Functions for repeat calcs

def CalcV(D):
    MFlowrate = (1.45 * 1E6) / 86400
    VFlowrate = MFlowrate/1200
    v = VFlowrate / (0.39 * (np.pi * D * D) / 4)
    return v

def Rate(Cnvrsn):
    RateA = 1+Cnvrsn**alpha
    RateB = (1-Cnvrsn)**beta
    return k*RateA*RateB

def CnvrsonArr(Bins):
    for m in range(Bins):
        Cnvrsn_m = Cnvrsn_Array[m]
        Cnvrsn_m = (dz/v)*Rate(Cnvrsn_m)+Cnvrsn_m
        if (Cnvrsn_m>1.):
            Cnvrsn_m=1.
        if (Cnvrsn_m<0.):
            Cnvrsn_m=0.
        Cnvrsn_Array[m+1] = Cnvrsn_m
    return Cnvrsn_Array

def Graph(z_Array, Cnvrsn_Array):
    plt.plot(z_Array,Cnvrsn_Array)
    plt.xlabel("z [m]")
    plt.ylabel("Cnvrsn")
    plt.axis([0.,np.amax(z_Array),0,1])
    plt.xticks(np.arange(np.amin(z_Array), np.amax(z_Array) + dz, 6*dz))
    plt.show()
    return

def Table(Nodes, z_Array, Cnvrsn_Array):
    print('D = {}m'.format(D))
    Line = "{:>12s}\t{:>12s}\t".format("z","Cnvrsn")
    print(Line)
    for m in range(Nodes):
        Line = "{:12.6f}\t{:12.6f}\t".format(z_Array[m],Cnvrsn_Array[m])
        print(Line)
    print('Surface area of reactor: {:.2f}m\n'.format(TSA))
    return

#Execution block

#D = 0.42 [m]
dz = 0.14
z_Array = np.arange(0,9.38+dz,dz)
Cnvrsn_Array = np.zeros_like(z_Array)
Nodes = len(z_Array)
Bins = Nodes-1
Cnvrsn_Array[0] = 0.01

D = 0.42 # [m]
v = CalcV(D)
Cnvrsn_Array = CnvrsonArr(Bins)
TSA = 2 *np.pi * D/2 * z_Array[-1]
Graph(z_Array, Cnvrsn_Array)
Table(Nodes, z_Array, Cnvrsn_Array)

#D = 0.84 [m]
dz = 0.14
z_Array = np.arange(0,2.24+dz,dz)
Cnvrsn_Array = np.zeros_like(z_Array)
Nodes = len(z_Array)
Bins = Nodes-1
Cnvrsn_Array[0] = 0.01

D = 0.84 # [m]
v = CalcV(D)
Cnvrsn_Array = CnvrsonArr(Bins)
TSA = 2 *np.pi * D/2 * z_Array[-1]
Graph(z_Array, Cnvrsn_Array)
Table(Nodes, z_Array, Cnvrsn_Array)

#D = 1.26 [m]
dz = 0.14/10
z_Array = np.arange(0,1.12+dz,dz)
Cnvrsn_Array = np.zeros_like(z_Array)
Nodes = len(z_Array)
Bins = Nodes-1
Cnvrsn_Array[0] = 0.01

D = 1.26 # [m]
v = CalcV(D)
Cnvrsn_Array = CnvrsonArr(Bins)
TSA = 2 *np.pi * D/2 * z_Array[-1]
Graph(z_Array, Cnvrsn_Array)
Table(Nodes, z_Array, Cnvrsn_Array)

#D = 1.68 [m]
dz = 0.14/10
z_Array = np.arange(0,0.7+dz,dz)
Cnvrsn_Array = np.zeros_like(z_Array)
Nodes = len(z_Array)
Bins = Nodes-1
Cnvrsn_Array[0] = 0.01

D = 1.68 # [m]
v = CalcV(D)
Cnvrsn_Array = CnvrsonArr(Bins)
TSA = 2 *np.pi * D/2 * z_Array[-1]
Graph(z_Array, Cnvrsn_Array)
Table(Nodes, z_Array, Cnvrsn_Array)




