# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 23:07:07 2021

@author: ianna
"""

import numpy as np
import matplotlib.pyplot as plt
#import time

k = 0.0135 # [/s] 
alpha = 0.24
beta = 1.87

VFlowrate = ((1.45 * 1E6) / 86400) / 1200

#Functions for repeat calcs

def CalcV(D):
    MFlowrate = (1.45 * 1E6) / 86400
    VFlowrate = MFlowrate/1200
    v = VFlowrate / (0.72 * (np.pi * D * D) / 4)
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
    print('Surface area of reactor: {:.2f}m\u00b2'.format(TSA))
    print('Volumetric flowrate of the polyepoxide: {:.4f}m\u00b3/s\n'.format(VFlowrate))
    return

#Execution block

#D = 0.42 [m]
dz = 0.14
z_Array = np.arange(0,76.58+dz,dz)
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
z_Array = np.arange(0,18.9+dz,dz)
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
dz = 0.14
z_Array = np.arange(0,8.26+dz,dz)
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
dz = 0.14
z_Array = np.arange(0,4.48+dz,dz)
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