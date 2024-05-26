# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 14:37:45 2022

@author: Alberto


second virial coefficient calculator = második viriális együttható kalkulátor

viriális for short.

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
filename=input('enter the name of the file:  ') 
#filename = 'h2-br2-cp-rydberg-cm.dat'  #butao.fname  #'H2F2.dat'#input('aaaaenter the name of the file:  ') 
arq = np.loadtxt(filename)



#CONSTANTS



rk = 7.24356  # 7.24297  #1/(8.6173324e-5)
h = 1.05459  # 1.05451   # 6.582119569e-16 #6.582119569  * 10e-16
N_A = 6.02252e-1
irk = 1.0 / rk





#LISTS CONTAINING ENTERED DATA

H = arq[0]
X = arq[1]
Z = arq[2]
Ta = arq[3]
Tb = arq[4]
L = arq[5]


H_a1 = float(H[0])
H_a2 = float(H[1])
H_a3 = float(H[2])
H_a4 = float(H[3])
H_a5 = float(H[4])
H_De = float(H[5])
H_Req = float(H[6])
H_Eref = float(H[7])

X_a1 = float(X[0])
X_a2 = float(X[1])
X_a3 = float(X[2])
X_a4 = float(X[3])
X_a5 = float(X[4])
X_De = float(X[5])
X_Req = float(X[6])
X_Eref = float(X[7])

Z_a1 = float(Z[0])
Z_a2 = float(Z[1])
Z_a3 = float(Z[2])
Z_a4 = float(Z[3])
Z_a5 = float(Z[4])
Z_De = float(Z[5])
Z_Req = float(Z[6])
Z_Eref = float(Z[7])

Ta_a1 = float(Ta[0])
Ta_a2 = float(Ta[1])
Ta_a3 = float(Ta[2])
Ta_a4 = float(Ta[3])
Ta_a5 = float(Ta[4])
Ta_De = float(Ta[5])
Ta_Req = float(Ta[6])
Ta_Eref = float(Ta[7])

Tb_a1 = float(Tb[0])
Tb_a2 = float(Tb[1])
Tb_a3 = float(Tb[2])
Tb_a4 = float(Tb[3])
Tb_a5 = float(Tb[4])
Tb_De = float(Tb[5])
Tb_Req = float(Tb[6])
Tb_Eref = float(Tb[7])

L_a1 = float(L[0])
L_a2 = float(L[1])
L_a3 = float(L[2])
L_a4 = float(L[3])
L_a5 = float(L[4])
L_De = float(L[5])
L_Req = float(L[6])
L_Eref = float(L[7])






#NUMERICAL DIFFERENTIATION
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




massaC = 12.0
massaO = 15.99491
massaN = 14.00307
massaH = 1.00783
massaF = 18.99840
massaBr = 79.904
massaCl = 34.96885
m1 = massaH
m2 = massaH
m3 = massaF #CL
m4 = massaF #Cl

#Br
B_ref = []
T_ref = []
#Cl
#B_ref = [-66.0103, -52.9951, -42.88,-34.849, -28.3552, -23.0219,  -18.5833,  -14.8471, -11.6713, -8.94904,  -6.59832, -4.55533,-2.76969, -1.20109,  0.18304, 1.40931,  2.49968,  3.47237,  4.34266, 5.12343, 5.82562, 6.45854,  7.0302, 7.54749, 8.01638,  8.44207, 8.82907, 9.18137, 9.50244,  9.79537,  10.0629,  10.3074,  10.531]
#T_ref = [300, 325, 350, 375,400, 425, 450, 475, 500, 525, 550, 575,600, 625, 650, 675,700, 725, 750, 775,800, 825, 850, 875, 900, 925, 950, 975,1000, 1025,  1050,  1075, 1100,]


#  H2F2
T_ref = [80, 90, 100, 110,120, 160, 200, 250, 300]
B_ref = [-50.13, -50.13, -33.91 ,-27.86 ,  -22.83 ,  -9.18,  -1.30,  4.70 , 8.49]




#MASS INPUT FOR MI

'''
m1 = input('qual atomo sera usado na primeira massa?')
m2 = input('qual atomo sera usado na segunda massa?')
m3 = input('qual atomo sera usado na terceira massa?')
m4 = input('qual atomo sera usado na quarta massa?')

'''
'''
B_ref = []
if m1 == 'O':
    m1 = 15.99491
    if m2 == 'O':
        T_ref = [70  ,85 ,100   ,120 ,140,170 ,200,240 ,280,330 ,380  ,430,495]
        B_ref = [-369.9,-262.9 ,-195.2,-138.8  ,-103.5  ,-69.7, -49.9, -32.2, -20.3, -10.1 , -2.9 , 2.6 ,8.0] 
    
elif m1 == 'F':
    m1 = 18.99840
    if m2 == 'F':
        T_ref = [85,95,105,120,135,150,170,190,210,230,250,270,295]
        B_ref = [-212.6, -172.1 ,-142.4 ,-110.3 ,-87.7 ,-71.0 ,-54.6 ,-42.5 ,-33.2 ,-25.9 ,-20.1 ,-15.2 ,-10.3]
        
elif m1 == 'C':
    m1 = 12.0
    if m2 == 'O':
        T_ref = [125,135,145,165,185,205,235,265,295,325,355,395,435,475,515,555,570]
        B_ref = [-118.1,-101.3,-87.6,-66.5,-51.2,-39.4,-26.3,-16.6,-9.1,-3.2,1.6,6.8,10.9,14.3,17.1,19.5,20.3]
elif m1 == 'N':
    m1 = 14.00307
    
    if m2 == 'O':
        T_ref = [125,135,150,165,185,215,245,285,325,375,425,475]
        B_ref = [-201.5, -158.1, -115.9, -89.4  , -67.1, -47.5, -35.2, -23.8, -15.1, -6.2, 1.3, 7.8]
    elif m2 == 'N':
        T_ref = [75,85,100,120,140,170,200,240,280,330,380,440,500,570,640,745]
        B_ref = [-276.8, -218.0, -160.7, -113.6 , -83.4, -54.4, -35.9, -19.6, -8.8, 0.5, 6.9, 12.4, 16.4, 19.8,22.5,25.3]
elif m1 == 'H':
    m1 = 1.00783
    if m2 == 'H':
       T2_ref = [15,20,25,30,35,40,45,50,55,60,65,75,85,100,115,130,145,165,180,200,220,250,290,330,370,410,450,490]
       T_ref = x = np.arange(60, 500, 10)
       B_ref = 1.7472*10 - (1.2926* 10**2)/(x) - (2.6988* 10**5)/(x)**2 + (8.0282 * 10**6)/(x)**3

       


       B2_ref = [-239.2,-150.6,-105.7,-79.0,-61.4,-49.0,-39.8,-32.7,-27.1,-22.6,-19.2,-13.2,-8.3,-2.8,1.2,4.2,6.4,8.6,9.8,11.1,12.1,13.2,14.1,14.8,15.3,15.7,15.9,16.2]

if m2=='O':
    m2=15.99491
elif m2=='F':
    m2= 18.99840
elif m2=='C':
    m2= 12.0
elif m2=='N':
    m2= 14.00307
elif m2=='H':
    m2= 1.00783
    
'''





Req1 = 0.744013
Req2 = 1.415962 # Cl  2.007880 # F = 1.415962 #usar dados CBS da tabela

mi = (m1 + m2) * (m3 + m4) / ((m1 + m2 + m3 + m4) * 6.02252)
mi1 = (m1 + m2)  /4 * Req1**2 / 6.02252
mi2 = (m3 + m4)  /4  * Req2**2 / 6.02252

print('mi1 mi2', mi1 , mi2)




    

# fit a straight line to the economic data

#B_ref = [-66.0103, -52.9951, -42.88,-34.849, -28.3552, -23.0219,  -18.5833,  -14.8471, -11.6713, -8.94904,  -6.59832, -4.55533,-2.76969, -1.20109,  0.18304, 1.40931,  2.49968,  3.47237,  4.34266, 5.12343, 5.82562, 6.45854,  7.0302, 7.54749, 8.01638,  8.44207, 8.82907, 9.18137, 9.50244,  9.79537,  10.0629,  10.3074,  10.531]
#T_ref = [300, 325, 350, 375,400, 425, 450, 475, 500, 525, 550, 575,600, 625, 650, 675,700, 725, 750, 775,800, 825, 850, 875, 900, 925, 950, 975,1000, 1025,  1050,  1075, 1100,]




def BH2_ref(T): # from 60K to 500K
    return 1.7472*10 - (1.2926* 10**2)/(T) - (2.6988* 10**5)/(T)**2 + (8.0282 * 10**6)/(T)**3
B_H2_ref = []
T_H2_ref = []   
for T in range(300,1000,5):
    B_H2_ref.append(BH2_ref(T))
    T_H2_ref.append(T)
    


def BCl2_ref(T): # from  300K  to 1070?
    return 1.3171*10 + 5.1055*10**3/T - 2.9404*10**7/T**2 - 5.0383*10**8/T**3

B_Cl2_ref = []
T_Cl2_ref = []

for T in range(300,1000,5):
    B_Cl2_ref.append(BCl2_ref(T))
    T_Cl2_ref.append(T)


def BF2_ref(T): # from  85K  to 295?
    return 3.3609*10 - 1.0625*10**4/T - 6.0780*10**5/T**2 - 2.2759*10**7/T**3

B_F2_ref = []
T_F2_ref = []

for T in range(300,1000,5):
    B_F2_ref.append(BF2_ref(T))
    T_F2_ref.append(T)
    
T_Br2_ref = [200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,950,1000]


B_Br2_ref = [-2161, -1869, -1637, -1449, -1294, -1165, -1056, -963.5, -883.4, -813.9, -753.0, -699.3, -651.6, -609.0, -570.8, -536.3, -505.1, -476.7, -450.7, -426.9, -405.0, -366.1, -332.7, -303.5, -278.0, -255.5, -235.4, -217.4, -201.2, -186.6, -173.3, -161.1, -150.0, -139.7, -130.2, -121.5, -113.4, -105.8, -98.7, -92.1, -85.9, -80.1, -74.6, -69.5, -64.6, -60.0, -49.5, -40.2]

# define the true objective function
'''
def objective(x, b1, b2, b3, b4, b5, b6):
    return b1 + b2/x + b3/x**2 + b4/x**3 + b5/x**4 + b6/x**5

# choose the input and output variables
# curve fit
popt1, _ = curve_fit(objective, T_H2_ref, B_H2_ref)
# summarize the parameter values
b1, b2, b3, b4, b5, b6 = popt1
print('y = %.3f + %.3fx^-1 + %.3fx^-2 + %.3fx^-3 + %.3fx^-4 + %.3fx^-5 ' % (b1, b2, b3, b4, b5, b6))
# plot input vs output
plt.scatter(T_H2_ref, B_H2_ref, label = 'ajuste H2')
# define a sequence of inputs between the smallest and largest known inputs
x_line = np.arange(min(T_H2_ref), max(T_H2_ref), 1)
# calculate the output for the range
y_line = objective(x_line, b1, b2, b3, b4, b5, b6)
# create a line plot for the mapping function
plt.plot(x_line, y_line, '--', color='red')

popt2, _ = curve_fit(objective, T_Cl2_ref, B_Cl2_ref)
# summarize the parameter values
b1, b2, b3, b4, b5, b6 = popt2
print('y = %.3f + %.3fx^-1 + %.3fx^-2 + %.3fx^-3 + %.3fx^-4 + %.3fx^-5 ' % (b1, b2, b3, b4, b5, b6))
# plot input vs output
plt.scatter(T_Cl2_ref, B_Cl2_ref, label = 'ajuste Cl2')
# define a sequence of inputs between the smallest and largest known inputs
x_line = np.arange(min(T_Cl2_ref), max(T_Cl2_ref), 1)
# calculate the output for the range
y_line = objective(x_line, b1, b2, b3, b4, b5, b6)
# create a line plot for the mapping function
plt.plot(x_line, y_line, '--', color='red')


popt3, _ = curve_fit(objective, T_F2_ref, B_F2_ref)
# summarize the parameter values
b1, b2, b3, b4, b5, b6 = popt3
print('y = %.3f + %.3fx^-1 + %.3fx^-2 + %.3fx^-3 + %.3fx^-4 + %.3fx^-5 ' % (b1, b2, b3, b4, b5, b6))
# plot input vs output
plt.scatter(T_F2_ref, B_F2_ref, label = 'ajuste F2')
# define a sequence of inputs between the smallest and largest known inputs
x_line = np.arange(min(T_F2_ref), max(T_F2_ref), 1)
# calculate the output for the range
y_line = objective(x_line, b1, b2, b3, b4, b5, b6)
# create a line plot for the mapping function
plt.plot(x_line, y_line, '--', color='red')
'''








'''
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

'''
#EQUAÇÃO VIRIAL DE ESTADO PARA MISTURAS 

Ma = m1 + m2 #verify
Mb = m3 + m4# esses dados sao do h2f2
w1 = -0.216
w2 = 0.053
rho_1 = 0.026667
rho_2 = 0.02715
zc1 = 0.303
zc2 = 0.287
x1 = 0.5
x2 = 0.5


#BH2_ref_list =[]
#BCl2_ref_list =[]

B_virial_state_ref = []
T_state = []
#
for T in range(300,1000,5):


   
    B_cruzado = (0.25*(22.1*Mb + Ma*(22.1 + Mb*T)) * (0.083 + 0.0695*w1 + 0.0695*w2) * (rho_1**(1/3) + rho_2**(1/3))**3) / ((10.9*Ma + 10.9*Mb + Ma*Mb*T)*(zc1 + zc2) * rho_1 * rho_2)

    B_virial_state = 2*x1*x2*B_cruzado + x2**2*BCl2_ref(T) + x1**2*BH2_ref(T)
    
    B_virial_state_ref.append(B_virial_state)
    T_state.append(T)




    plt.plot(T_state, B_virial_state_ref, '--', color='black')


plt.legend()
plt.show()


def UH(r):     
        y=(r-H_Req)/H_Req
        UH = -H_De*(1 + H_a1*y + H_a2*y**2 + H_a3*y**3 + H_a4*y**4 + H_a5*y**5) * np.exp(-H_a1*y) + H_Eref
        
      #esse potencial está certo
        #print(r, UH)
        return UH 
     
def UX(r):

        y=(r-X_Req)/X_Req
        UX = -X_De*(1 + X_a1*y + X_a2*y**2 + X_a3*y**3 + X_a4*y**4 + X_a5*y**5) * np.exp(-X_a1*y) + X_Eref

        
      #esse potencial está certo
 
       
        #print(r, UX)
        return UX
    
def UZ(r):

        y=(r-Z_Req)/Z_Req
        UZ = -Z_De*(1 + Z_a1*y + Z_a2*y**2 + Z_a3*y**3 + Z_a4*y**4 + Z_a5*y**5) * np.exp(-Z_a1*y) + Z_Eref

      #esse potencial está certo
 
       
        #print(r, UZ)
        return UZ
    
def UTa(r):
        y=(r-Ta_Req)/Ta_Req
        UTa = -Ta_De*(1 + Ta_a1*y + Ta_a2*y**2 + Ta_a3*y**3 + Ta_a4*y**4 + Ta_a5*y**5) * np.exp(-Ta_a1*y) + Ta_Eref

       
        #print(r, UTa)
        return UTa
    
def UTb(r):
        y=(r-Tb_Req)/Tb_Req
        UTb = -Tb_De*(1 + Tb_a1*y + Tb_a2*y**2 + Tb_a3*y**3 + Tb_a4*y**4 + Tb_a5*y**5) * np.exp(-Tb_a1*y) + Tb_Eref

      #esse potencial está certo
 
       
        #print(r, UTb)
        return UTb
 
def UL(r):

        y=(r-L_Req)/L_Req
        UL = -L_De*(1 + L_a1*y + L_a2*y**2 + L_a3*y**3 + L_a4*y**4 + L_a5*y**5) * np.exp(-L_a1*y) + L_Eref

      #esse potencial está certo
 
       
        #print(r, UL)
        return UL


def UM_000(r):
    UM_000=(2*UH(r)+UL(r)+2*(UTa(r)+UTb(r)+UX(r)))/9
    return UM_000

def UM_202(r):
    UM_202=-2*(UH(r)-UL(r)+UTa(r)-2*UTb(r)+UX(r))/(9*(5**(1/2)))    
    return UM_202
def UM_022(r):
    UM_022=-2*(UH(r)-UL(r)-2*UTa(r)+UTb(r)+UX(r))/(9*(5**(1/2)))    
    return UM_022
def UM_220(r):
    UM_220= 2*(4*UH(r)-UL(r)-5*UTa(r)-5*UTb(r)-5*UX(r)+12*UZ(r))/(45*(5**(1/2)))
    return UM_220
def UM_222(r):
    UM_222=((2/7)**(1/2))*(13*UH(r)-UL(r)+7*UTa(r)+7*UTb(r)-14*UX(r)-12*UZ(r))/45  
    return UM_222

def UM_224(r):
    UM_224=((2/35)**(1/2))*8*(UH(r)+UL(r)-2*UZ(r))/15    
    return UM_224


def UM_FINAL(r, th_a, th_b, phi): #potencial que vai ser usado na equação 3
    UM_FINAL = 1e-3*(UM_000(r)+(5**(1/2))/4 * (3*(np.cos(2*th_a))+1) * UM_202(r) 
    + (5**(1/2))/4 * (3*(np.cos(2*th_b))+1) * UM_022(r) 
    + (5**(1/2))/16 * UM_220(r) * ((3*(np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
    + 12*np.sin(2*th_a)* np.sin(2*th_b)*np.cos(phi)
    + 3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)) 
    - 14**(1/2)*5/112 * UM_222(r) * ((3*np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
    + 6*np.sin(2*th_a)*np.sin(2*th_b) * np.cos(phi) 
    - 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b))*np.cos(2*phi)
    + (3*(70)**(1/2))/112 * UM_224(r) * ((3*(np.cos(2*th_a))+1)*(3*(np.cos(2*th_b))+1) 
    - 8* np.sin(2*th_a) * np.sin(2*th_b)*np.cos(phi)+((1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi))/2))
    #print('r:',r ,'000:',UM_000(r),'202:',UM_202(r),'022:',UM_022(r),'220:',UM_220(r),'222:',UM_222(r),'224:',UM_224(r))
    #print('um',UM_FINAL)
    return UM_FINAL


def d1_H(r):
    d1_H=derivative(UH,r,method='central',p=1e-8)
    return d1_H
def d1_L(r):
    d1_L=derivative(UL,r,method='central',p=1e-8)
    return d1_L
def d1_Ta(r): 
    d1_Ta=derivative(UTa,r,method='central',p=1e-8)
    return d1_Ta
def d1_Tb(r):
    d1_Tb=derivative(UTb,r,method='central',p=1e-8)
    return d1_Tb
def d1_X(r):
    d1_X=derivative(UX,r,method='central',p=1e-8)
    return d1_X
def d1_Z(r):
    d1_Z=derivative(UZ,r,method='central',p=1e-8)
    return d1_Z
def d2_H(r):
    d2_H=derivative2(UH,r,method='central',p=1e-8)
    return d2_H
def d2_L(r):
    d2_L=derivative2(UL,r,method='central',p=1e-8)
    return d2_L
def d2_Ta(r):
    d2_Ta=derivative2(UTa,r,method='central',p=1e-8)
    return d2_Ta
def d2_Tb(r):
    d2_Tb=derivative2(UTb,r,method='central',p=1e-8)
    return d2_Tb
def d2_X(r):
    
    d2_X=derivative2(UX,r,method='central',p=1e-8)
    return d2_X
def d2_Z(r):
    d2_Z=derivative2(UZ,r,method='central',p=1e-8)
    return d2_Z

def d1r(r, th_a, th_b, phi):
    d1r = (-1/18)*(1+3*np.cos(2*th_a))*(d1_H(r)-d1_L(r)+d1_Ta(r)-2*d1_Tb(r)+d1_X(r))
    -(1/18)*(1+3*np.cos(2*th_b))*(d1_H(r)-d1_L(r)-2*d1_Ta(r)+d1_Tb(r)+d1_X(r))
    +(1/9)*(2*d1_H(r)+d1_L(r)+2*(d1_Ta(r)+d1_Tb(r)+d1_X(r)))
    -(1/504)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))
    -3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)
    +6*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(13*d1_H(r)
    -d1_L(r)+7*d1_Ta(r)+7*d1_Tb(r)-14*d1_X(r)-12*d1_Z(r))
    +1/35*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+1/2*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)
    -8*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(d1_H(r)+d1_L(r)-2*d1_Z(r))
    +1/360*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)
    +12*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(4*d1_H(r)-d1_L(r)-5*d1_Ta(r)-5*d1_Tb(r)-5*d1_X(r)+12*d1_Z(r))
    return d1r

def d1th_a(r, th_a, th_b, phi):
    d1th_a = 1/3*np.sin(2*th_a) * (UH(r)-UL(r)+UTa(r)-2*UTb(r)+UX(r))
    -1/504*(-6*(1+3*np.cos(2*th_b)) * np.sin(2*th_a) -6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
    +12*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (13*UH(r)-UL(r)+7*UTa(r)+7*UTb(r)-14*UX(r)-12*UZ(r))
    +1/35*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + (1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
    -16*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (UH(r)+UL(r)-2*UZ(r))
    +1/360*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + 6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
    +24*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (4*UH(r)-UL(r)-5*UTa(r)-5*UTb(r)-5*UX(r)+12*UZ(r))
    
    return d1th_a


def d1th_b(r, th_a, th_b, phi):
    d1th_b = 1/3*np.sin(2*th_b) * (UH(r)-UL(r)-2*UTa(r)+UTb(r)+UX(r))
    -1/504*(12*np.cos(2*th_b)* np.cos(phi) * np.sin(2*th_a) -6*(1+3*np.cos(2*th_a)) * np.sin(2*th_b) -6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(13*UH(r)-UL(r)+7*UTa(r)+7*UTb(r)-14*UX(r)-12*UZ(r))
    +1/35*(-16*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a)
    -6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b) + (1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(UH(r)+UL(r)-2*UZ(r))
    +1/360*(24*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a) -6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b)
    + 6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))* (4*UH(r)-UL(r)-5*UTa(r)-5*UTb(r)-5*UX(r)+12*UZ(r))
    return d1th_b


def d1phi(r, th_a, th_b, phi):
    d1phi = -1/504*(-6*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) + 6*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi))* (13*UH(r)-UL(r)+7*UTa(r)+7*UTb(r)-14*UX(r)-12*UZ(r))
    +1/35*(8*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - (1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi))* (UH(r)+UL(r)-2*UZ(r))
    + 1/360*(-12*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - 6*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi))*(4*UH(r)-UL(r)-5*UTa(r)-5*UTb(r)-5*UX(r)+12*UZ(r))
    return d1phi


def d2r(r, th_a, th_b, phi):
    d2r = -1/18*(1+3*np.cos(2*th_a)) * (d2_H(r)-d2_L(r)+d2_Ta(r)-2*d2_Tb(r)+d2_X(r))
    -1/18*(1+3*np.cos(2*th_b)) * (d2_H(r)-d2_L(r)-2*d2_Ta(r)+d2_Tb(r)+d2_X(r))
    +1/9* (2*d2_H(r)+d2_L(r)+2*(d2_Ta(r)+d2_Tb(r)+d2_X(r)))
    -1/504*(1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b))
    -3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi) + 6*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b) * (13*d2_H(r)-d2_L(r)+7*d2_Ta(r)+7*d2_Tb(r)-14*d2_X(r)-12*d2_Z(r))
    +1/35*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 1/2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
    -8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (d2_H(r)+d2_L(r)-2*d2_Z(r))
    +1/360*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
    +12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*d2_H(r)-d2_L(r)-5*d2_Ta(r)-5*d2_Tb(r)-5*d2_X(r)+12*d2_Z(r))
    return d2r


def d2th_a(r, th_a, th_b, phi):
    d2th_a = 2/3*np.cos(2*th_a) * (UH(r)-UL(r)+UTa(r)-2*UTb(r)+UX(r))
    + 1/504*(12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
    + 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*UH(r)-UL(r)+7*UTa(r)+7*UTb(r)-14*UX(r)-12*UZ(r))
    + 1/35*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 2*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
    + 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (UH(r)+UL(r)-2*UZ(r))
    + 1/360*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
    - 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*UH(r)-UL(r)-5*UTa(r)-5*UTb(r)-5*UX(r)+12*UZ(r))
    return d2th_a


def d2th_b(r, th_a, th_b, phi):
    d2th_b = 2/3*np.cos(2*th_b) * (UH(r)-UL(r)-2*UTa(r)+UTb(r)+UX(r))
    + 1/504*(12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 12*(1-np.cos(2*th_a)) * np.cos(2*th_b) * np.cos(2*phi)
    + 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*UH(r)-UL(r)+7*UTa(r)+7*UTb(r)-14*UX(r)-12*UZ(r))
    + 1/35*(-12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 2*(1-np.cos(2*th_a) * np.cos(2*th_b)) * np.cos(2*phi)
    + 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (UH(r)+UL(r)-2*UZ(r))
    + 1/360*(-12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 12*(1-np.cos(2*th_a)) * np.cos(2*th_b) * np.cos(2*phi)
    - 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*UH(r)-UL(r)-5*UTa(r)-5*UTb(r)-5*UX(r)+12*UZ(r))
    return d2th_b


def d2phi(r, th_a, th_b, phi):
    d2phi = 1/504*(-12*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
    + 6*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*UH(r)-UL(r)+7*UTa(r)+7*UTb(r)-14*UX(r)-12*UZ(r))
    + 1/35*(-2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
    + 8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (UH(r)+UL(r)-2*UZ(r))
    + 1/360*(-12*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
    - 12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*UH(r)-UL(r)-5*UTa(r)-5*UTb(r)-5*UX(r)+12*UZ(r))
    return d2phi

Reqs = np.array([H_Req, X_Req, Z_Req, Ta_Req, Tb_Req, L_Req])


lims = mean(Reqs)
lim_inf = lims/4

lim_sup = 10*lims/2
print(lim_inf, lim_sup)


r0 = 4.22 #mean(Reqs)

def integrand_vegas(x):
    r    = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*rk/T
    if F > 700:
        F = 700
    #print('f',F,x)
    

    return  N_A * np.pi * np.sin(th_a) * r0**3/3 * np.sin(th_b)*(r**2)*(1-np.exp(-F))#*10e-24 #24 eh o certo de ang^3 pra cm ^3




def integrand_c1(x): #primeira correção quantica
    r     = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*rk/T
    if F > 700:
        F = 700
    d_1r= d1r(r, th_a, th_b, phi)
    #c1 = N_A * h**2 * rk**3 / (48 * mi * T**3) * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * d_1r**2
    c1 = np.pi/24 * N_A * (h**2)/mi *  (rk/T)**3  * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * d_1r**2 * r**2

    

    
    return c1

def integrand_c2(x): #segunda correção quantica 
    r     = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*rk/T
    if F > 700:
        F = 700
    d_1th_a = d1th_a(r, th_a, th_b, phi)
    d_1th_b = d1th_b(r, th_a, th_b, phi)
    d_2th_a = d2th_a(r, th_a, th_b, phi)
    d_2th_b = d2th_b(r, th_a, th_b, phi)
    d_2phi = d2phi(r, th_a, th_b, phi)
    
    
    
    
    VL1 = -d_1th_a * np.cos(th_a)/np.sin(th_a) - d_2th_a - d_2phi/(np.sin(th_a))**2
    VL2 = -d_1th_b * np.cos(th_b)/np.sin(th_b) - d_2th_b - d_2phi/(np.sin(th_b))**2
    
    c2 = -np.pi/24 * N_A * h**2 * r**2 * np.exp(-F) * (VL1/mi1 + VL2/mi2)  * (rk/T)**2
    
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
    F=UM_FINAL(r,th_a,th_b,phi)*rk/T
    if F > 700:
        F = 700
    d_1th_a = d1th_a(r, th_a, th_b, phi)
    d_1th_b = d1th_b(r, th_a, th_b, phi)
    d_2th_a = d2th_a(r, th_a, th_b, phi)
    d_2th_b = d2th_b(r, th_a, th_b, phi)
    d_2phi = d2phi(r, th_a, th_b, phi)
    d_1r= d1r(r, th_a, th_b, phi)
    d_2r = d2r(r, th_a, th_b, phi)
    
    #c3a =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/(np.power(np.sin(th_a),2)) * d_2phi) * r**2
    #c3b =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/(np.power(np.sin(th_b),2)) * d_2phi) * r**2
    
    c3a_max = np.exp(-F) * r**2 * (d_2r**2/10 + d_1r**2/(5*r**2) + d_1r**3/(9/rk*r*T) - d_1r**4/(72*(1/rk*T)**2))
    #c3b_max = 
    
    #c3a = np.sin(th_a) * np.sin(th_b) * np.exp(-F)
    #c3b = (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/(np.power(np.sin(th_a),2)) * d_2phi)
    #c3c = (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/(np.power(np.sin(th_b),2)) * d_2phi) * r**2
    c3 = -np.pi/48 * N_A * (h**4)/mi**2 * (c3a_max ) * (rk/T)**4 * np.sin(th_a) * np.sin(th_b)
    #talves ese rk / T ou nao #*( h**2)/(2 * mi * 1**2)
    


    return c3



def integrand_c4(x): #quarta correção quantica
    r     = x[0]
    th_a = x[1]
    th_b = x[2]
    phi  = x[3] 
    F=UM_FINAL(r,th_a,th_b,phi)*rk/T
    if F > 700:
        F = 700 
    d_1r = d1r(r, th_a, th_b, phi)
    d_2r = d2r(r, th_a, th_b, phi)
    d_1th_a = d1th_a(r, th_a, th_b, phi)
    d_1th_b = d1th_b(r, th_a, th_b, phi)
    d_2th_a = d2th_a(r, th_a, th_b, phi)
    d_2th_b = d2th_b(r, th_a, th_b, phi)
    d_2phi = d2phi(r, th_a, th_b, phi)
    

    VL1 = -d_1th_a * np.cos(th_a)/np.sin(th_a) - d_2th_a - d_2phi/(np.sin(th_a))**2
    VL2 = -d_1th_b * np.cos(th_b)/np.sin(th_b) - d_2th_b - d_2phi/(np.sin(th_b))**2
    
    
    c4 = (VL1 + VL2) * np.exp(-F) * (-np.pi/24) * N_A * h**2/mi * (rk/T)**2



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


#minha integração monte carlo
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
    B_clas.append(result.mean ) 
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

#f, (ax1, ax2,ax3, ax4) = plt.subplots(1, 4, constrained_layout=True)

f, (ax1, ax2,ax3) = plt.subplots(1, 3, constrained_layout=True)
ax1.plot(r, UH(r), label = 'UH')
ax1.plot(r, UX(r), label = 'UX')
ax1.plot(r, UZ(r), label = 'UZ')
ax1.plot(r, UTa(r), label = 'UTa')
ax1.plot(r, UTb(r), label = 'UTb')
ax1.plot(r, UL(r), label = 'UL')
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



ax2.plot(r, UM_000(r),  color='black',label = 'U000')
ax2.plot(r, UM_202(r), color='green',label = 'U202')
ax2.plot(r, UM_220(r), color='blue',label = 'U220')
ax2.plot(r, UM_022(r), color='red',label = 'U022')
ax2.plot(r, UM_222(r), color='cyan',label = 'U222')
ax2.plot(r, UM_224(r), color='magenta',label = 'U224')
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

constantes = []
aa = 0
for y in Tstr:
    for i in T_ref:
        if i == y:
            index_i = T_ref.index(i)
            index_y = Tstr.index(y)
            aa+=1    
            #print(Bsem[index_y], y, index_y)
            #print(B_ref[index_i], i, index_i)
            print('B_ref=', B_ref[index_i],'Bmeu=', B_clas[index_y],'T=', y)
            
            cte = B_clas[index_y]/B_ref[index_i]
            
            constantes.append(cte)
            
print(constantes)





i = 1
def objective(T, a1, a2, a3, a4, a5 , a6):
    
    ajuste = a1  + a2  * T**(-1) + a3 * T**(-2) + a4 * T**(-3) + a5 * T**(-4)+ a6 * T**(-5) 
    
    return ajuste
    #função a_i x^(-i+1)






popt, _ = curve_fit(objective, Tstr, B_plus_all_except_c2)

a1, a2, a3, a4, a5, a6 = popt

print(a1, a2, a3, a4, a5, a6) 
ax3.plot(Tstr,B_plus_all_except_c2, '.', label = 'pontos')

x_line = arange(min(Tstr), max(Tstr), 1)

y_line = objective(x_line, a1, a2, a3, a4, a5, a6)


ax3.plot(x_line, y_line, color = 'pink',label = 'Bfit')
ax3.legend()
pyplot.show()



export_list = [[B_plus_all_except_c2], [Tstr]]





df = pd.DataFrame(export_list)
writer = pd.ExcelWriter('h2brz2export.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name='welcome', index=False)
writer.save()


#print('r:',r ,'000:',UM_000(r),'202:',UM_202(r),'022:',UM_022(r),'220:',UM_220(r),'222:',UM_222(r),'224:',UM_224(r))
