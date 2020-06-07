# -*- coding: utf-8 -*-
"""
Created on Mon May  4 17:02:24 2020

@author: vondreele
this is a python translation of a fortran routine in RFMCProfile that does the
convolution of F(Q) and sinc(Q*Rmax) where Rmax is radius of largest sphere in
big box model. The result is input into fullrmc as ReducedStructureFactor
"""
import numpy as np
def sincConvolution(XY,d):

    n = XY.shape[1]
    snew = np.zeros(n)
    dq = np.zeros(n)
    sold = XY[1]
    q = XY[0]
    dq[1:] = np.diff(q)
    dq[0] = dq[1]
    
    for j in range(n):
        for i in range(n):
            b = abs(q[i]-q[j])
            t = q[i]+q[j]
            if j == i:
                snew[j] += q[i]*sold[i]*(d-np.sin(t*d)/t)*dq[i]
            else:
                snew[j] += q[i]*sold[i]*(np.sin(b*d)/b-np.sin(t*d)/t)*dq[i]
        snew[j] /= np.pi*q[j]
    
    snew[0] = snew[1]
    return snew
