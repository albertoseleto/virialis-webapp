# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 16:44:29 2022

@author: Alberto
"""

import vegas
import random
import time
import math
import os
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean
import pandas as pd
from numpy import arange
from pandas import read_csv
from matplotlib import pyplot



#TIME MEASURE

start = time.time()


#DATA ENTRY
#filename=input('enter the name of the file:  ') 

data = pd.read_csv('H2F2.dat', sep="\s+", header=None)
data.columns = ["a1", "a2", "a3", "a4", 'a5', 'De', 'Req', 'Eref', 'err1', 'err2']

print(data) #Dtype = float64, acho que não vai dar problema



#CONSTANTS

class constants:

    rk = 7.24356  # 7.24297  #1/(8.6173324e-5)
    h = 1.05459  # 1.05451   # 6.582119569e-16 #6.582119569  * 10e-16
    N_A = 6.02252e-1
    irk = 1.0 / rk

#NUMERICAL DIFFERENTIATION

class num_diff:
    
    def derivative(g,a,method='central',p=0.01):
        '''Compute the difference formula for f'(a) with step size h.
    
        Parameters
        ----------
        f : function
            Vectorized function of one variable
        a : number
            Compute derivative at x = a
        method : string
            Difference formula: 'forward', 'backward' or 'central'
        h : number
            Step size in difference formula
    
        Returns
        -------
        float
            Difference formula:
                central: f(a+h) - f(a-h))/2h
                forward: f(a+h) - f(a))/h
                backward: f(a) - f(a-h))/h            
        '''
        if method == 'central':
            return (g(a + p) - g(a - p))/(2*p)
        elif method == 'forward':
            return (g(a + p) - g(a))/p
        elif method == 'backward':
            return (g(a) - g(a - p))/p
        else:
            raise ValueError("Method must be 'central', 'forward' or 'backward'.")
                 
    def derivative2(g,a,method='central',p=0.01): #derivada de segunda ordem
    
    
        '''Compute the difference formula for f'(a) with step size h.
    
        Parameters
        ----------
        f : function
            Vectorized function of one variable
        a : number
            Compute derivative at x = a
        method : string
            Difference formula: 'forward', 'backward' or 'central'
        h : number
            Step size in difference formula
    
        Returns
        -------
        float
            Difference formula:
                central: f(a+h) - f(a-h))/2h
                forward: f(a+h) - f(a))/h
                backward: f(a) - f(a-h))/h            
        '''
        if method == 'central':
            return (g(a + p) - 2*g(a)+g(a - p))/(p**2)
        
        elif method == 'forward':
            return (g(a + 2*p) -2*g(a+p)+ g(a))/p
        
        elif method == 'backward':
            return (g(a) -2*g(a - p)+g(a-2*p))/p
        else:
            raise ValueError("Method must be 'central', 'forward' or 'backward'.")

#CLASSES CONTAINING ENTERED DATA

class H: 
    
    a1 = data.loc[0,'a1']
    a2 = data.loc[0,'a2']
    a3 = data.loc[0,'a3']
    a4 = data.loc[0,'a4']
    a5 = data.loc[0,'a5']   
    De = data.loc[0,'De']
    Req = data.loc[0,'Req']
    Eref = data.loc[0,'Eref']
      
class X: 
    
    a1 = data.loc[1,'a1']
    a2 = data.loc[1,'a2']
    a3 = data.loc[1,'a3']
    a4 = data.loc[1,'a4']
    a5 = data.loc[1,'a5']   
    De = data.loc[1,'De']
    Req = data.loc[1,'Req']
    Eref = data.loc[1,'Eref']
       
class Z: 
    
    a1 = data.loc[2,'a1']
    a2 = data.loc[2,'a2']
    a3 = data.loc[2,'a3']
    a4 = data.loc[2,'a4']
    a5 = data.loc[2,'a5']   
    De = data.loc[2,'De']
    Req = data.loc[2,'Req']
    Eref = data.loc[2,'Eref']
    
class Ta: 
    
    a1 = data.loc[3,'a1']
    a2 = data.loc[3,'a2']
    a3 = data.loc[3,'a3']
    a4 = data.loc[3,'a4']
    a5 = data.loc[3,'a5']   
    De = data.loc[3,'De']
    Req = data.loc[3,'Req']
    Eref = data.loc[3,'Eref']
    
class Tb: 
    
    a1 = data.loc[4,'a1']
    a2 = data.loc[4,'a2']
    a3 = data.loc[4,'a3']
    a4 = data.loc[4,'a4']
    a5 = data.loc[4,'a5']   
    De = data.loc[4,'De']
    Req = data.loc[4,'Req']
    Eref = data.loc[4,'Eref']
    

class L: 
    
    a1 = data.loc[5,'a1']
    a2 = data.loc[5,'a2']
    a3 = data.loc[5,'a3']
    a4 = data.loc[5,'a4']
    a5 = data.loc[5,'a5']   
    De = data.loc[5,'De']
    Req = data.loc[5,'Req']
    Eref = data.loc[5,'Eref']
    
    
class LC:
    def UH(r):     
        y=(r-H.Req)/H.Req
        UH = -H.De*(1 + H.a1*y + H.a2*y**2 + H.a3*y**3 + H.a4*y**4 + H.a5*y**5) * np.exp(-H.a1*y) + H.Eref
        return UH 
     
    def UX(r):
    
        y=(r-X.Req)/X.Req
        UX = -X.De*(1 + X.a1*y + X.a2*y**2 + X.a3*y**3 + X.a4*y**4 + X.a5*y**5) * np.exp(-X.a1*y) + X.Eref 
        return UX
        
    def UZ(r):
    
        y=(r-Z.Req)/Z.Req
        UZ = -Z.De*(1 + Z.a1*y + Z.a2*y**2 + Z.a3*y**3 + Z.a4*y**4 + Z.a5*y**5) * np.exp(-Z.a1*y) + Z.Eref        
        return UZ
        
    def UTa(r):
        y=(r-Ta.Req)/Ta.Req
        UTa = -Ta.De*(1 + Ta.a1*y + Ta.a2*y**2 + Ta.a3*y**3 + Ta.a4*y**4 + Ta.a5*y**5) * np.exp(-Ta.a1*y) + Ta.Eref
        return UTa
        
    def UTb(r):
        y=(r-Tb.Req)/Tb.Req
        UTb = -Tb.De*(1 + Tb.a1*y + Tb.a2*y**2 + Tb.a3*y**3 + Tb.a4*y**4 + Tb.a5*y**5) * np.exp(-Tb.a1*y) + Tb.Eref
        return UTb
     
    def UL(r):
        y=(r-L.Req)/L.Req
        UL = -L.De*(1 + L.a1*y + L.a2*y**2 + L.a3*y**3 + L.a4*y**4 + L.a5*y**5) * np.exp(-L.a1*y) + L.Eref
        return UL
    
 
class complexa:
    def UM_000(r):
        UM_000=(2*LC.UH(r)+LC.UL(r)+2*(LC.UTa(r)+LC.UTb(r)+LC.UX(r)))/9
        return UM_000
    
    def UM_202(r):
        UM_202=-2*(LC.UH(r)-LC.UL(r)+LC.UTa(r)-2*LC.UTb(r)+LC.UX(r))/(9*(5**(1/2)))    
        return UM_202
    def UM_022(r):
        UM_022=-2*(LC.UH(r)-LC.UL(r)-2*LC.UTa(r)+LC.UTb(r)+LC.UX(r))/(9*(5**(1/2)))    
        return UM_022
    def UM_220(r):
        UM_220= 2*(4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))/(45*(5**(1/2)))
        return UM_220
    def UM_222(r):
        UM_222=((2/7)**(1/2))*(13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))/45  
        return UM_222
    
    def UM_224(r):
        UM_224=((2/35)**(1/2))*8*(LC.UH(r)+LC.UL(r)-2*LC.UZ(r))/15    
        return UM_224
         
def UM_FINAL(r, th_a, th_b, phi): #potencial que vai ser usado na equação 3
    UM_FINAL = 1e-3*(complexa.UM_000(r)+(5**(1/2))/4 * (3*(np.cos(2*th_a))+1) * complexa.UM_202(r) 
    + (5**(1/2))/4 * (3*(np.cos(2*th_b))+1) * complexa.UM_022(r) 
    + (5**(1/2))/16 * complexa.UM_220(r) * ((3*(np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
    + 12*np.sin(2*th_a)* np.sin(2*th_b)*np.cos(phi)
    + 3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)) 
    - 14**(1/2)*5/112 * complexa.UM_222(r) * ((3*np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
    + 6*np.sin(2*th_a)*np.sin(2*th_b) * np.cos(phi) 
    - 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b))*np.cos(2*phi)
    + (3*(70)**(1/2))/112 * complexa.UM_224(r) * ((3*(np.cos(2*th_a))+1)*(3*(np.cos(2*th_b))+1) 
    - 8* np.sin(2*th_a) * np.sin(2*th_b)*np.cos(phi)+((1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi))/2))
    #print('r:',r ,'000:',UM_000(r),'202:',UM_202(r),'022:',UM_022(r),'220:',UM_220(r),'222:',UM_222(r),'224:',UM_224(r))
    #print('um',UM_FINAL)
    return UM_FINAL




class drv:
        
    def d1_H(r):
        d1_H = num_diff.derivative(LC.UH,r,method='central',p=1e-8)
        return d1_H
    def d1_L(r):
        d1_L = num_diff.derivative(LC.UL,r,method='central',p=1e-8)
        return d1_L
    def d1_Ta(r): 
        d1_Ta = num_diff.derivative(LC.UTa,r,method='central',p=1e-8)
        return d1_Ta
    def d1_Tb(r):
        d1_Tb = num_diff.derivative(LC.UTb,r,method='central',p=1e-8)
        return d1_Tb
    def d1_X(r):
        d1_X = num_diff.derivative(LC.UX,r,method='central',p=1e-8)
        return d1_X
    def d1_Z(r):
        d1_Z = num_diff.derivative(LC.UZ,r,method='central',p=1e-8)
        return d1_Z
    def d2_H(r):
        d2_H = num_diff.derivative2(LC.UH,r,method='central',p=1e-8)
        return d2_H
    def d2_L(r):
        d2_L = num_diff.derivative2(LC.UL,r,method='central',p=1e-8)
        return d2_L
    def d2_Ta(r):
        d2_Ta = num_diff.derivative2(LC.UTa,r,method='central',p=1e-8)
        return d2_Ta
    def d2_Tb(r):
        d2_Tb = num_diff.derivative2(LC.UTb,r,method='central',p=1e-8)
        return d2_Tb
    def d2_X(r):
        d2_X = num_diff.derivative2(LC.UX,r,method='central',p=1e-8)
        return d2_X
    def d2_Z(r):
        d2_Z = num_diff.derivative2(LC.UZ,r,method='central',p=1e-8)
        return d2_Z
    
    def d1r(r, th_a, th_b, phi):
        d1r = (-1/18)*(1+3*np.cos(2*th_a))*(drv.d1_H(r)-drv.d1_L(r)+drv.d1_Ta(r)-2*drv.d1_Tb(r)+drv.d1_X(r))
        -(1/18)*(1+3*np.cos(2*th_b))*(drv.d1_H(r)-drv.d1_L(r)-2*drv.d1_Ta(r)+drv.d1_Tb(r)+drv.d1_X(r))
        +(1/9)*(2*drv.d1_H(r)+drv.d1_L(r)+2*(drv.d1_Ta(r)+drv.d1_Tb(r)+drv.d1_X(r)))
        -(1/504)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))
        -3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)
        +6*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(13*drv.d1_H(r)
        -drv.d1_L(r)+7*drv.d1_Ta(r)+7*drv.d1_Tb(r)-14*drv.d1_X(r)-12*drv.d1_Z(r))
        +1/35*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+1/2*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)
        -8*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(drv.d1_H(r)+drv.d1_L(r)-2*drv.d1_Z(r))
        +1/360*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)
        +12*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(4*drv.d1_H(r)-drv.d1_L(r)-5*drv.d1_Ta(r)-5*drv.d1_Tb(r)-5*drv.d1_X(r)+12*drv.d1_Z(r))
        return d1r
    
    def d1th_a(r, th_a, th_b, phi):
        d1th_a = 1/3*np.sin(2*th_a) * (LC.UH(r)-LC.UL(r)+LC.UTa(r)-2*LC.UTb(r)+LC.UX(r))
        -1/504*(-6*(1+3*np.cos(2*th_b)) * np.sin(2*th_a) -6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
        +12*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
        +1/35*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + (1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
        -16*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
        +1/360*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + 6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
        +24*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
        
        return d1th_a
    
    
    def d1th_b(r, th_a, th_b, phi):
        d1th_b = 1/3*np.sin(2*th_b) * (LC.UH(r)-LC.UL(r)-2*LC.UTa(r)+LC.UTb(r)+LC.UX(r))
        -1/504*(12*np.cos(2*th_b)* np.cos(phi) * np.sin(2*th_a) -6*(1+3*np.cos(2*th_a)) * np.sin(2*th_b) -6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
        +1/35*(-16*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a)
        -6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b) + (1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
        +1/360*(24*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a) -6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b)
        + 6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))* (4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
        return d1th_b
    
    
    def d1phi(r, th_a, th_b, phi):
        d1phi = -1/504*(-6*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) + 6*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi))* (13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
        +1/35*(8*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - (1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi))* (LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
        + 1/360*(-12*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - 6*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi))*(4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
        return d1phi
    
    
    def d2r(r, th_a, th_b, phi):
        d2r = -1/18*(1+3*np.cos(2*th_a)) * (drv.d2_H(r)-drv.d2_L(r)+drv.d2_Ta(r)-2*drv.d2_Tb(r)+drv.d2_X(r))
        -1/18*(1+3*np.cos(2*th_b)) * (drv.d2_H(r)-drv.d2_L(r)-2*drv.d2_Ta(r)+drv.d2_Tb(r)+drv.d2_X(r))
        +1/9* (2*drv.d2_H(r)+drv.d2_L(r)+2*(drv.d2_Ta(r)+drv.d2_Tb(r)+drv.d2_X(r)))
        -1/504*(1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b))
        -3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi) + 6*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b) * (13*drv.d2_H(r)-drv.d2_L(r)+7*drv.d2_Ta(r)+7*drv.d2_Tb(r)-14*drv.d2_X(r)-12*drv.d2_Z(r))
        +1/35*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 1/2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
        -8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (drv.d2_H(r)+drv.d2_L(r)-2*drv.d2_Z(r))
        +1/360*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
        +12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*drv.d2_H(r)-drv.d2_L(r)-5*drv.d2_Ta(r)-5*drv.d2_Tb(r)-5*drv.d2_X(r)+12*drv.d2_Z(r))
        return d2r
    
    
    def d2th_a(r, th_a, th_b, phi):
        d2th_a = 2/3*np.cos(2*th_a) * (LC.UH(r)-LC.UL(r)+LC.UTa(r)-2*LC.UTb(r)+LC.UX(r))
        + 1/504*(12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
        + 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
        + 1/35*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 2*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
        + 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
        + 1/360*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
        - 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
        return d2th_a
    
    
    def d2th_b(r, th_a, th_b, phi):
        d2th_b = 2/3*np.cos(2*th_b) * (LC.UH(r)-LC.UL(r)-2*LC.UTa(r)+LC.UTb(r)+LC.UX(r))
        + 1/504*(12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 12*(1-np.cos(2*th_a)) * np.cos(2*th_b) * np.cos(2*phi)
        + 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
        + 1/35*(-12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 2*(1-np.cos(2*th_a) * np.cos(2*th_b)) * np.cos(2*phi)
        + 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
        + 1/360*(-12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 12*(1-np.cos(2*th_a)) * np.cos(2*th_b) * np.cos(2*phi)
        - 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
        return d2th_b
    
    
    def d2phi(r, th_a, th_b, phi):
        d2phi = 1/504*(-12*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
        + 6*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
        + 1/35*(-2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
        + 8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
        + 1/360*(-12*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
        - 12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
        return d2phi
 

class B_references:
    
    def BH2_ref(T): # from 60K to 500K
        return 1.7472*10 - (1.2926* 10**2)/(T) - (2.6988* 10**5)/(T)**2 + (8.0282 * 10**6)/(T)**3
    
    def BCl2_ref(T): # from  300K  to 1070?
        return 1.3171*10 + 5.1055*10**3/T - 2.9404*10**7/T**2 - 5.0383*10**8/T**3

    def BF2_ref(T): # from  85K  to 295?
        return 3.3609*10 - 1.0625*10**4/T - 6.0780*10**5/T**2 - 2.2759*10**7/T**3
    
    
massaC = 12.0
massaO = 15.99491
massaN = 14.00307
massaH = 1.00783
massaF = 18.99840
massaBr = 79.904
massaCl = 34.96885            
      

B_A2_ref = []
T_A2_ref = []   

     
B_B2_ref = []
T_B2_ref = []


      
while True:
    molecula = str(input('qual molécula será estudada? (formato de input: A2B2):'))
    if molecula == 'H2F2':
        m1 = massaH
        m2 = massaH
        m3 = massaF 
        m4 = massaF 
        
        Req1 = 0.744013
        Req2 = 1.415962 # Cl  2.007880 # F = 1.415962 #usar dados CBS da tabela
        T_ref = [80, 90, 100, 110,120, 160, 200, 250, 300]
        B_ref = [-50.13, -50.13, -33.91 ,-27.86 ,  -22.83 ,  -9.18,  -1.30,  4.70 , 8.49] 
        
        
        mi = (m1 + m2) * (m3 + m4) / ((m1 + m2 + m3 + m4) * 6.02252)
        mi1 = (m1 + m2)  /4 * Req1**2 / 6.02252
        mi2 = (m3 + m4)  /4  * Req2**2 / 6.02252

        Ma = m1 + m2 
        Mb = m3 + m4
        w1 = -0.216
        w2 = 0.053
        rho_1 = 0.026667
        rho_2 = 0.02715
        zc1 = 0.303
        zc2 = 0.287
        x1 = 0.5
        x2 = 0.5
        def B_A2_ref(T):
            return B_references.BF2_ref(T)
        def B_B2_ref(T):
            return B_references.BH2_ref(T)
        break
        
    elif molecula == 'H2H2':
        m1 = massaH
        m2 = massaH
        m3 = massaH
        m4 = massaH     
        Req1 = 0.744013
        Req2 = 0.744013
        
        mi = (m1 + m2) * (m3 + m4) / ((m1 + m2 + m3 + m4) * 6.02252)
        mi1 = (m1 + m2)  /4 * Req1**2 / 6.02252
        mi2 = (m3 + m4)  /4  * Req2**2 / 6.02252
        
        Ma = m1 + m2 
        Mb = m3 + m4
        w1 = -0.216
        w2 = -0.216
        rho_1 = 0.026667
        rho_2 = 0.026667
        zc1 = 0.303
        zc2 = 0.303
        x1 = 0.5
        x2 = 0.5
        
        def B_A2_ref(T):
            return B_references.BH2_ref(T)
        def B_B2_ref(T):
            return B_references.BH2_ref(T)
        
        
        break
    
    elif molecula == 'H2Cl2':
        m1 = massaH
        m2 = massaH
        m3 = massaCl
        m4 = massaCl  
        B_ref = [-66.0103, -52.9951, -42.88,-34.849, -28.3552, -23.0219,  -18.5833,  -14.8471, -11.6713, -8.94904,  -6.59832, -4.55533,-2.76969, -1.20109,  0.18304, 1.40931,  2.49968,  3.47237,  4.34266, 5.12343, 5.82562, 6.45854,  7.0302, 7.54749, 8.01638,  8.44207, 8.82907, 9.18137, 9.50244,  9.79537,  10.0629,  10.3074,  10.531]
        T_ref = [300, 325, 350, 375,400, 425, 450, 475, 500, 525, 550, 575,600, 625, 650, 675,700, 725, 750, 775,800, 825, 850, 875, 900, 925, 950, 975,1000, 1025,  1050,  1075, 1100,]
        Req1 = 0.744013
        Req2 = 2.007880 
        
        mi = (m1 + m2) * (m3 + m4) / ((m1 + m2 + m3 + m4) * 6.02252)
        mi1 = (m1 + m2)  /4 * Req1**2 / 6.02252
        mi2 = (m3 + m4)  /4  * Req2**2 / 6.02252
        
        Ma = m1 + m2 #verify
        Mb = m3 + m4# esses dados sao do h2cl2
        w1 = -0.216
        w2 = 0.279
        rho_1 = 0.026667
        rho_2 = 0.05087
        zc1 = 0.303
        zc2 = 0.073
        x1 = 0.5
        x2 = 0.5
        
        def B_A2_ref(T):
            return B_references.BCl2_ref(T)
        def B_B2_ref(T):
            return B_references.BH2_ref(T)
        
        

        break
    elif molecula == 'H2Br2':
        m1 = massaH
        m2 = massaH
        m3 = massaBr
        m4 = massaBr
        B_ref = []
        T_ref = []
        Req1 = 0
        Req2 = 0
        
        
        mi = (m1 + m2) * (m3 + m4) / ((m1 + m2 + m3 + m4) * 6.02252)
        mi1 = (m1 + m2)  /4 * Req1**2 / 6.02252
        mi2 = (m3 + m4)  /4  * Req2**2 / 6.02252
        
        Ma = m1 + m2 #verify
        Mb = m3 + m4# esses dados sao do h2br2
        w1 = -0.216
        w2 = 0.286
        rho_1 = 0.026667
        rho_2 = 0.05538
        zc1 = 0.303
        zc2 = 0.126
        x1 = 0.5
        x2 = 0.5

        break
        
    else: 
        print('atualmente, temos suporte as moléculas H2F2, H2H2, e H2Cl2. Tente uma delas')
        continue
            

B_virial_state_ref = []
T_state = []


for T in range(300,1000,5):


   
    B_cruzado = (0.25*(22.1*Mb + Ma*(22.1 + Mb*T)) * (0.083 + 0.0695*w1 + 0.0695*w2) * (rho_1**(1/3) + rho_2**(1/3))**3) / ((10.9*Ma + 10.9*Mb + Ma*Mb*T)*(zc1 + zc2) * rho_1 * rho_2)

    B_virial_state = 2*x1*x2*B_cruzado + x2**2*B_A2_ref(T) + x1**2*B_B2_ref(T)
    
    B_virial_state_ref.append(B_virial_state)
    T_state.append(T)

#plt.plot(T_state, B_virial_state_ref, '--', color='black')
#plt.legend()
#plt.show()


Reqs = np.array([H.Req, X.Req, Z.Req, Ta.Req, Tb.Req, L.Req])

lims = mean(Reqs)
lim_inf = lims/4
lim_sup = 10*lims/2
r0 = 4.22 #mean(Reqs)

def integrand_vegas(x):
    r    = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*constants.rk/T
    if F > 700:
        F = 700
    #print('f',F,x)
    

    return  constants.N_A * np.pi * np.sin(th_a) * r0**3/3 * np.sin(th_b)*(r**2)*(1-np.exp(-F))#*10e-24 #24 eh o certo de ang^3 pra cm ^3

def integrand_c1(x): #primeira correção quantica
    r     = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*constants.rk/T
    if F > 700:
        F = 700
    d_1r= drv.d1r(r, th_a, th_b, phi)
    #c1 = N_A * h**2 * rk**3 / (48 * mi * T**3) * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * d_1r**2
    c1 = np.pi/24 * constants.N_A * (constants.h**2)/mi *  (constants.rk/T)**3  * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * d_1r**2 * r**2

    return c1

def integrand_c2(x): #segunda correção quantica 
    r     = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*constants.rk/T
    if F > 700:
        F = 700
    d_1th_a = drv.d1th_a(r, th_a, th_b, phi)
    d_1th_b = drv.d1th_b(r, th_a, th_b, phi)
    d_2th_a = drv.d2th_a(r, th_a, th_b, phi)
    d_2th_b = drv.d2th_b(r, th_a, th_b, phi)
    d_2phi = drv.d2phi(r, th_a, th_b, phi)
    
    
    
    
    VL1 = -d_1th_a * np.cos(th_a)/np.sin(th_a) - d_2th_a - d_2phi/(np.sin(th_a))**2
    VL2 = -d_1th_b * np.cos(th_b)/np.sin(th_b) - d_2th_b - d_2phi/(np.sin(th_b))**2
    
    c2 = -np.pi/24 * constants.N_A * constants.h**2 * r**2 * np.exp(-F) * (VL1/mi1 + VL2/mi2)  * (constants.rk/T)**2
    
    #c2a =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/(np.power(np.sin(th_a),2)) * d_2phi) * r**2  #talvez sem esse rk/T
#-np.pi/12 * N_A * h**2/mi1 *
    #c2b =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/(np.power(np.sin(th_b),2)) * d_2phi) * r**2  #talvez sem esse rk/T

    #c2 = -np.pi/12 * N_A * h**2 * (c2a/mi1 +  c2b/mi2) * (rk/T)**2
    return c2

def integrand_c3(x): #terceira correção quantica nossa, segunda do maximiliano 
    r     = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*constants.rk/T
    if F > 700:
        F = 700
    d_1th_a = drv.d1th_a(r, th_a, th_b, phi)
    d_1th_b = drv.d1th_b(r, th_a, th_b, phi)
    d_2th_a = drv.d2th_a(r, th_a, th_b, phi)
    d_2th_b = drv.d2th_b(r, th_a, th_b, phi)
    d_2phi = drv.d2phi(r, th_a, th_b, phi)
    d_1r= drv.d1r(r, th_a, th_b, phi)
    d_2r = drv.d2r(r, th_a, th_b, phi)
    
    #c3a =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/(np.power(np.sin(th_a),2)) * d_2phi) * r**2
    #c3b =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/(np.power(np.sin(th_b),2)) * d_2phi) * r**2
    
    c3a_max = np.exp(-F) * r**2 * (d_2r**2/10 + d_1r**2/(5*r**2) + d_1r**3/(9/constants.rk*r*T) - d_1r**4/(72*(1/constants.rk*T)**2))
    #c3b_max = 
    
    #c3a = np.sin(th_a) * np.sin(th_b) * np.exp(-F)
    #c3b = (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/(np.power(np.sin(th_a),2)) * d_2phi)
    #c3c = (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/(np.power(np.sin(th_b),2)) * d_2phi) * r**2
    c3 = -np.pi/48 * constants.N_A * (constants.h**4)/mi**2 * (c3a_max ) * (constants.rk/T)**4 * np.sin(th_a) * np.sin(th_b)
    #talves ese rk / T ou nao #*( h**2)/(2 * mi * 1**2)
    return c3

def integrand_c4(x): #quarta correção quantica
    r     = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*constants.rk/T
    if F > 700:
        F = 700 
    d_1r = drv.d1r(r, th_a, th_b, phi)
    d_2r = drv.d2r(r, th_a, th_b, phi)
    d_1th_a = drv.d1th_a(r, th_a, th_b, phi)
    d_1th_b = drv.d1th_b(r, th_a, th_b, phi)
    d_2th_a = drv.d2th_a(r, th_a, th_b, phi)
    d_2th_b = drv.d2th_b(r, th_a, th_b, phi)
    d_2phi = drv.d2phi(r, th_a, th_b, phi)
    

    VL1 = -d_1th_a * np.cos(th_a)/np.sin(th_a) - d_2th_a - d_2phi/(np.sin(th_a))**2
    VL2 = -d_1th_b * np.cos(th_b)/np.sin(th_b) - d_2th_b - d_2phi/(np.sin(th_b))**2
    
    
    c4 = (VL1 + VL2) * np.exp(-F) * (-np.pi/24) * constants.N_A * constants.h**2/mi * (constants.rk/T)**2

    #c4a = -N_A * h**4 * rk**4 / (1920 * mi**2 * T**4) * np.sin(th_a) * np.sin(th_b) * np.exp(-F)
    #c4b = (d_2r**2 + 2/r**2*(d_1r**2) + 10 * rk/(9 * T * r) * d_1r**3 - 5 * rk**2/(36 * T**2) * d_1r**4) * r**2
    #c4aa = c4a * c4b
    #c4p = np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (d_2r**2 + 2/r**2*(d_1r**2) + 10 * rk/(9 * T * r) * d_1r**3 - 5 * rk**2/(36 * T**2) * d_1r**4) * r**2  #talvez sem rk/T
    #c4 = -np.pi/12 * N_A * h**2/mi * c4p *(rk/T)**2
    return c4


integ = vegas.Integrator([[lim_inf, lim_sup], [0, np.pi], [0, np.pi], [0, 2*np.pi]])

B_clas = []
Tstr = []

B_plus_c1=[]
B_main = []
B_c1 = []
B_c2 = []
B_c3 = []
B_c4 = []
B_correcoes = []
B_plus_all_except_c2 = []


for T in range(50,1000,100):   #

    integ(integrand_vegas, nitn=10, neval=1000)  

    integ(integrand_c1, nitn=10, neval=1000)  

    integ(integrand_c2, nitn=10, neval=1000)   
    integ(integrand_c3, nitn=10, neval=1000)   
    integ(integrand_c4, nitn=10, neval=1000)   

    result = integ(integrand_vegas, nitn=10, neval=1000)   
    
    
    print(result.summary())
    print('result VEGAS = %s    Q = %.2f' % (result, result.Q)) 
    print(T)
    B_clas.append(result.mean) 
    Tstr.append(T)
    
    result_c1 = integ(integrand_c1, nitn=10, neval=1000)   

    result_c2 = integ(integrand_c2, nitn=10, neval=1000)   
    result_c3 = integ(integrand_c3, nitn=10, neval=1000)   
    result_c4 = integ(integrand_c4, nitn=10, neval=1000)

    

    print('result VEGAS c1 = %s    Q = %.2f' % (result_c1, result_c1.Q)) 
    print('result VEGAS c2 = %s    Q = %.2f' % (result_c2, result_c2.Q)) 
    print('result VEGAS c3 = %s    Q = %.2f' % (result_c3, result_c3.Q)) 
    print('result VEGAS c4 = %s    Q = %.2f' % (result_c4, result_c4.Q)) 

    print('result VEGAS geral =' , result.mean + result_c1.mean + result_c2.mean + result_c3.mean + result_c4.mean)


    B_main.append(result.mean + result_c1.mean - result_c2.mean + result_c3.mean + result_c4.mean)

    

    B_c1.append(-result_c1.mean)
    B_c2.append(result_c2.mean)
    B_c3.append(result_c3.mean)
    B_c4.append(-result_c4.mean)
    B_correcoes.append(result_c2.mean + result_c3.mean + result_c1.mean + result_c4.mean)
    
    
    
    B_plus_all_except_c2.append(result.mean + result_c1.mean + result_c3.mean + result_c4.mean)
    
    #B_cruzado = (0.25*(22.1*Mb + Ma*(22.1 + Mb*T)) * (0.083 + 0.0695*w1 + 0.0695*w2) * (rho_1**(1/3) + rho_2**(1/3))**3) / ((10.9*Ma + 10.9*Mb + Ma*Mb*T)*(zc1 + zc2) * rho_1 * rho_2)

    #B_virial_state = 2*x1*x2*B_cruzado + x2**2*BCl2_ref(T) + x1**2*BH2_ref(T)
    
    #B_virial_state_ref.append(B_virial_state)
    
#graficos
r = np.linspace(0, 10, 100)
th_a = np.linspace(0, np.pi, 100)
th_b = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)


f, (ax1, ax2,ax3) = plt.subplots(1, 3, constrained_layout=True)
ax1.plot(r, LC.UH(r), label = 'UH')
ax1.plot(r, LC.UX(r), label = 'UX')
ax1.plot(r, LC.UZ(r), label = 'UZ')
ax1.plot(r, LC.UTa(r), label = 'UTa')
ax1.plot(r, LC.UTb(r), label = 'UTb')
ax1.plot(r, LC.UL(r), label = 'UL')
ax1.set_title('(a)energias simples')
#ax1.set_xlim([2, 8]) #1.2 B 
#ax1.set_ylim([-40, 100]) #inferior : B(T=10) * 1.1 superior: B(T=500)* 1.2
ax1.set_ylabel(r'$U(r)[eV]$')
ax1.set_xlabel(r'$r[Å]$', labelpad=1)
#ax2.set_xlim([2, 8]) #1.2 B 
#ax2.set_ylim([-40, 100])
#ax3.set_xlim([2, 10]) #1.2 B 
#ax3.set_ylim([-1000, 1006])

#ax4.set_xlim([2, 10]) #1.2 B 
#ax4.set_ylim([-1000, 1006])



ax2.plot(r, complexa.UM_000(r),  color='black',label = 'U000')
ax2.plot(r, complexa.UM_202(r), color='green',label = 'U202')
ax2.plot(r, complexa.UM_220(r), color='blue',label = 'U220')
ax2.plot(r, complexa.UM_022(r), color='red',label = 'U022')
ax2.plot(r, complexa.UM_222(r), color='cyan',label = 'U222')
ax2.plot(r, complexa.UM_224(r), color='magenta',label = 'U224')
ax2.set_title('(b)energias complexas')
#ax2.set_xlim([2, 9]) #1.2 B 
#ax2.set_ylim([-0.4, 2]) #inferior : B(T=10) * 1.1 superior: B(T=500)* 1.2
ax2.set_ylabel(r'$U(r)[eV]$')
ax2.set_xlabel(r'$r[Å]$', labelpad=1)

ax3.plot(Tstr, B_clas, color='r',label = 'Bclássico')
#ax3.plot(T_ref, B_ref,   color='b',label = 'ref') #Breferencia
ax3.plot(T_state, B_virial_state_ref,   color='g',label = 'Bvirialstate')



#ax3.plot(Tstr, B_main,   color='yellow',label = 'Bmain')
ax3.plot(Tstr, B_plus_all_except_c2, color='blue',label = 'Btotal')

ax3.set_title('(c) Segundo Coeficiente Virial(B) em função da Temperatura(T)')
ax3.set_ylabel(r'$B[cm^3/mol]$')
ax3.set_xlabel(r'$T[K]$', labelpad=1)
'''
ax4.plot(Tstr, B_c1, color='g',label = 'Bc1')
ax4.plot(Tstr, B_c2,  color='red',label = 'Bc2')
ax4.plot(Tstr, B_c3,   color='blue',label = 'Bc3')
ax4.plot(Tstr, B_c4,   color='magenta',label = 'Bc4')

#ax3.set_ylim([-100000, 100000]) #1.2 B 



ax4.plot(Tstr, B_correcoes,   color='black',label = 'Bcorreções')



ax_image = fig.add_subplot(Tstr, Bstr)
ax_image.set_title('Imagem original')
ax_image.imshow(image, cmap='gray')
'''

#plt.subplot(Tstr, Bstr)
#plt.scatter(Tstr, Bstr)
#plt.title(f'Gráficos (a) do , name. You are age.')

ax1.legend()
ax2.legend()
ax3.legend()
#ax4.legend()


plt.show()
print('evaluation time {}'.format(time.strftime("%H:%M:%S", time.gmtime(time.time()-start))))
print('number of temperatures:', len(Tstr))
