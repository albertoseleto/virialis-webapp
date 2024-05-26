# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:00:27 2023

@author: Alberto
"""

import streamlit as st
import numpy as np
import pandas as pd
import vegas
import time
from scipy.optimize import curve_fit
from statistics import mean
from numpy import arange
from pandas import read_csv
from plotly.subplots import make_subplots
import plotly.graph_objects as go





# Specify what pages should be shown in the sidebar, and what their titles 
# and icons should be



# TIME MEASURE





# DATA ENTRY
st.title('Second Virial Coefficient Calculator for A2B2 Molecule')

uploaded_file = st.file_uploader("upload a file")

potential = st.selectbox(
    'What potential energy do you want to use?',
    ('Rydberg Potential', 'Improved Leonard-Jonnes Potential'))


step = st.text_input(
    'For what Temperature do you want to calculate?')


try:
    step = float(step)
except:
    pass


gas_type = st.selectbox(
    'What gas structure are you studying?',
    ('A2B2', 'AB',))


st.write('The calculus will be made using ', potential,
         'on a gas of structure ', gas_type, 'using a temperature of', step)

if gas_type == 'A2B2':

    monomero1 = st.selectbox(
        'What is your first monomer?',
        ('H2', 'F2', 'Cl2', 'Br2'))
    monomero2 = st.selectbox(
        'What is your second monomer?',
        ('H2', 'F2', 'Cl2', 'Br2'))
    gas = monomero1 + monomero2

data_virial = pd.DataFrame()

class constants:

    rk = 7.24356  # 7.24297  #1/(8.6173324e-5)
    h = 1.05459  # 1.05451   # 6.582119569e-16 #6.582119569  * 10e-16
    N_A = 6.02252e-1
    irk = 1.0 / rk


class num_diff:

    def derivative(g, a, method='central', p=0.01):
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
            raise ValueError(
                "Method must be 'central', 'forward' or 'backward'.")

    def derivative2(g, a, method='central', p=0.01):  # derivada de segunda ordem
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
            return (g(a + 2*p) - 2*g(a+p) + g(a))/p

        elif method == 'backward':
            return (g(a) - 2*g(a - p)+g(a-2*p))/p
        else:
            raise ValueError(
                "Method must be 'central', 'forward' or 'backward'.")


if st.button('Calculate'):
    start = time.time()
    if uploaded_file is not None:

        if gas_type == 'A2B2':

            st.write('Calculating for ', gas, '...')

            if potential == 'Improved Leonard-Jonnes Potential' :
                st.write('LEONARD-JONNES IS IN CONSTRUCTION...')

                data = pd.read_csv(uploaded_file, sep="\s+", header=None)
                data.columns = ["alpha", "beta", "mp", "De", "Req" ]
                st.subheader('DataFrame')
                st.write(data)

                class H:

                    alpha = data.loc[0, 'alpha']
                    beta = data.loc[0, 'beta']
                    mp = data.loc[0, 'mp']
                    De = data.loc[0, 'De'] #* 0.000123981 #entrada em cm-1, passando para EV********VALOR CORRETO VERIFICADO **********(0.002859/23.06055)#* 0.0012398 #(0.002859/23.06055) / 1.2398e-7
                    Req = data.loc[0, 'Req']
                    

                class X:

                    alpha = data.loc[1, 'alpha']
                    beta = data.loc[1, 'beta']
                    mp = data.loc[1, 'mp']
                    De = data.loc[1, 'De'] #* 0.000123981 #entrada em cm-1, passando para EV********VALOR CORRETO VERIFICADO **********(0.002859/23.06055)#* 0.0012398 #(0.002859/23.06055) / 1.2398e-7
                    Req = data.loc[1, 'Req']
                    

                class Z:

                    alpha = data.loc[2, 'alpha']
                    beta = data.loc[2, 'beta']
                    mp = data.loc[2, 'mp']
                    De = data.loc[2, 'De'] #* 0.000123981 #entrada em cm-1, passando para EV********VALOR CORRETO VERIFICADO **********(0.002859/23.06055)#* 0.0012398 #(0.002859/23.06055) / 1.2398e-7
                    Req = data.loc[2, 'Req']
                    

                class Ta:

                    alpha = data.loc[3, 'alpha']
                    beta = data.loc[3, 'beta']
                    mp = data.loc[3, 'mp']
                    De = data.loc[3, 'De'] #* 0.000123981 #entrada em cm-1, passando para EV********VALOR CORRETO VERIFICADO **********(0.002859/23.06055)#* 0.0012398 #(0.002859/23.06055) / 1.2398e-7
                    Req = data.loc[3, 'Req']
                


                class Tb:

                    alpha = data.loc[4, 'alpha']
                    beta = data.loc[4, 'beta']
                    mp = data.loc[4, 'mp']
                    De = data.loc[4, 'De'] #* 0.000123981 #entrada em cm-1, passando para EV********VALOR CORRETO VERIFICADO **********(0.002859/23.06055)#* 0.0012398 #(0.002859/23.06055) / 1.2398e-7
                    Req = data.loc[4, 'Req']
                    


                class L:

                    alpha = data.loc[5, 'alpha']
                    beta = data.loc[5, 'beta']
                    mp = data.loc[5, 'mp']
                    De = data.loc[5, 'De'] #* 0.000123981 #entrada em cm-1, passando para EV********VALOR CORRETO VERIFICADO **********(0.002859/23.06055)#* 0.0012398 #(0.002859/23.06055) / 1.2398e-7
                    Req = data.loc[5, 'Req']
                    

                '''
                class Sa:

                    alpha = data.loc[6, 'alpha']
                    beta = data.loc[6, 'beta']
                    mp = data.loc[6, 'mp']
                    Req = data.loc[6, 'Req']
                    De = data.loc[6, 'De']  #* 0.000123981 #entrada em cmcm-1, passando para EV********VALOR CORRETO VERIFICADO **********(0.002859/23.06055)#* 0.0012398 #(0.002859/23.06055) / 1.2398e-7



                class Sb:

                    alpha = data.loc[7, 'alpha']
                    beta = data.loc[7, 'beta']
                    mp = data.loc[7, 'mp']
                    Req = data.loc[7, 'Req']
                    De = data.loc[7, 'De'] #* 0.000123981 #entrada em cm-1, passando para EV********VALOR CORRETO VERIFICADO **********(0.002859/23.06055)#* 0.0012398 #(0.002859/23.06055) / 1.2398e-7
                '''


                class LC:  # formula ta certa
                    def UH(r):

                        n = H.beta + H.alpha * (r / H.Req) ** 2
                        return H.De * ((H.mp/(n - H.mp) * (H.Req/r) ** n) - (n/(n - H.mp) * (H.Req/r) ** H.mp))

                    def UX(r):

                        n = X.beta + X.alpha * (r / X.Req) ** 2
                        return X.De * ((X.mp/(n - X.mp) * (X.Req/r) ** n) - (n/(n - X.mp) * (X.Req/r) ** X.mp))

                    def UZ(r):

                        n = Z.beta + Z.alpha * (r / Z.Req) ** 2
                        return Z.De * ((Z.mp/(n - Z.mp) * (Z.Req/r) ** n) - (n/(n - Z.mp) * (Z.Req/r) ** Z.mp))

                    def UTa(r):
                        n = Ta.beta + Ta.alpha * (r / Ta.Req) ** 2
                        return Ta.De * ((Ta.mp/(n - Ta.mp) * (Ta.Req/r) ** n) - (n/(n - Ta.mp) * (Ta.Req/r) ** Ta.mp))

                    def UTb(r):
                        n = Tb.beta + Tb.alpha * (r / Tb.Req) ** 2
                        return Tb.De * ((Tb.mp/(n - Tb.mp) * (Tb.Req/r) ** n) - (n/(n - Tb.mp) * (Tb.Req/r) ** Tb.mp))

                    def UL(r):
                        n = L.beta + L.alpha * (r / L.Req) ** 2
                        return L.De * ((L.mp/(n - L.mp) * (L.Req/r) ** n) - (n/(n - L.mp) * (L.Req/r) ** L.mp))
                    '''
                    def USa(r):
                        n = Sa.beta + Sa.alpha * (r / Sa.Req) ** 2
                        return Sa.De * ((Sa.mp/(n - Sa.mp) * (Sa.Req/r) ** n) - (n/(n - Sa.mp) * (Sa.Req/r) ** Sa.mp))

                    def USb(r):
                        n = Sb.beta + Sb.alpha * (r / Sb.Req) ** 2
                        return Sb.De * ((Sb.mp/(n - Sb.mp) * (Sb.Req/r) ** n) - (n/(n - Sb.mp) * (Sb.Req/r) ** Sb.mp))
                    '''

                class complexa:
                    def UM_000(r):
                        UM_000 = (2*LC.UH(r)+LC.UL(r) + 2 *
                                  (LC.UTa(r)+LC.UTb(r)+LC.UX(r)))/9
                        return UM_000

                    def UM_202(r):
                        UM_202 = 2*(LC.UH(r)-LC.UL(r)+LC.UTa(r) -
                                    2*LC.UTb(r)+LC.UX(r))/(9*(5**(1/2)))
                        return UM_202

                    def UM_022(r):
                        UM_022 = 2*(LC.UH(r)-LC.UL(r)-2*LC.UTa(r) +
                                    LC.UTb(r)+LC.UX(r))/(9*(5**(1/2)))
                        return UM_022

                    def UM_220(r):
                        UM_220 = 2*(4*LC.UH(r)-LC.UL(r)-5*(LC.UTa(r) + 
                                    LC.UTb(r) + LC.UX(r))+12*LC.UZ(r))/(45*(5**(1/2)))
                        return UM_220

                    def UM_222(r):
                        UM_222 = ((2/7)**(1/2))*(13*LC.UH(r)-LC.UL(r)+7 *
                                                 (LC.UTa(r)+LC.UTb(r)-2*LC.UX(r))-12*LC.UZ(r))/45
                        return UM_222

                    def UM_224(r):
                        UM_224 = ((2/35)**(1/2)*8 * LC.UH(r)+LC.UL(r)+2*LC.UZ(r))/15
                        return UM_224

                def UM_FINAL(r, th_a, th_b, phi):  # potencial que vai ser usado na equação 3
                    UM_FINAL = (complexa.UM_000(r)+(5**(1/2))/4 * (3*(np.cos(2*th_a))+1) * complexa.UM_202(r) 
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
                    # print('um',UM_FINAL)
                    return UM_FINAL

            elif potential == 'Rydberg Potential':

                data = pd.read_csv(uploaded_file, sep="\s+", header=None)

                data.columns = ["a1", "a2", "a3", "a4", 'a5',
                                'De', 'Req', 'Eref']
                st.subheader('DataFrame')
                st.write(data)
                #st.subheader('Descriptive Statistics')
                # st.write(data.describe())

                class H:

                    a1 = data.loc[0, 'a1']
                    a2 = data.loc[0, 'a2']
                    a3 = data.loc[0, 'a3']
                    a4 = data.loc[0, 'a4']
                    a5 = data.loc[0, 'a5'] 
                    De = data.loc[0, 'De'] 
                    Req = data.loc[0, 'Req']
                    Eref = data.loc[0, 'Eref']

                class X:

                    a1 = data.loc[1, 'a1']
                    a2 = data.loc[1, 'a2']
                    a3 = data.loc[1, 'a3']
                    a4 = data.loc[1, 'a4']
                    a5 = data.loc[1, 'a5']
                    De = data.loc[1, 'De']  # eV
                    Req = data.loc[1, 'Req']
                    Eref = data.loc[1, 'Eref']

                class Z:

                    a1 = data.loc[2, 'a1']
                    a2 = data.loc[2, 'a2']
                    a3 = data.loc[2, 'a3']
                    a4 = data.loc[2, 'a4']
                    a5 = data.loc[2, 'a5']
                    De = data.loc[2, 'De']
                    Req = data.loc[2, 'Req']
                    Eref = data.loc[2, 'Eref']

                class Ta:

                    a1 = data.loc[3, 'a1']
                    a2 = data.loc[3, 'a2']
                    a3 = data.loc[3, 'a3']
                    a4 = data.loc[3, 'a4']
                    a5 = data.loc[3, 'a5']
                    De = data.loc[3, 'De']
                    Req = data.loc[3, 'Req']
                    Eref = data.loc[3, 'Eref']

                class Tb:

                    a1 = data.loc[4, 'a1']
                    a2 = data.loc[4, 'a2']
                    a3 = data.loc[4, 'a3']
                    a4 = data.loc[4, 'a4']
                    a5 = data.loc[4, 'a5']
                    De = data.loc[4, 'De']
                    Req = data.loc[4, 'Req']
                    Eref = data.loc[4, 'Eref']

                class L:

                    a1 = data.loc[5, 'a1']
                    a2 = data.loc[5, 'a2']
                    a3 = data.loc[5, 'a3']
                    a4 = data.loc[5, 'a4']
                    a5 = data.loc[5, 'a5']
                    De = data.loc[5, 'De']
                    Req = data.loc[5, 'Req']
                    Eref = data.loc[5, 'Eref']

                class LC:
                    def UH(r):
                        y = (r-H.Req)/H.Req
                        UH = -H.De*(1 + H.a1*y + H.a2*y**2 + H.a3*y**3 +
                                    H.a4*y**4 + H.a5*y**5) * np.exp(-H.a1*y) + H.Eref
                        return UH

                    def UX(r):

                        y = (r-X.Req)/X.Req
                        UX = -X.De*(1 + X.a1*y + X.a2*y**2 + X.a3*y**3 +
                                    X.a4*y**4 + X.a5*y**5) * np.exp(-X.a1*y) + X.Eref
                        return UX

                    def UZ(r):

                        y = (r-Z.Req)/Z.Req
                        UZ = -Z.De*(1 + Z.a1*y + Z.a2*y**2 + Z.a3*y**3 +
                                    Z.a4*y**4 + Z.a5*y**5) * np.exp(-Z.a1*y) + Z.Eref
                        return UZ

                    def UTa(r):
                        y = (r-Ta.Req)/Ta.Req
                        UTa = -Ta.De*(1 + Ta.a1*y + Ta.a2*y**2 + Ta.a3*y**3 +
                                      Ta.a4*y**4 + Ta.a5*y**5) * np.exp(-Ta.a1*y) + Ta.Eref
                        return UTa

                    def UTb(r):
                        y = (r-Tb.Req)/Tb.Req
                        UTb = -Tb.De*(1 + Tb.a1*y + Tb.a2*y**2 + Tb.a3*y**3 +
                                      Tb.a4*y**4 + Tb.a5*y**5) * np.exp(-Tb.a1*y) + Tb.Eref
                        return UTb

                    def UL(r):
                        y = (r-L.Req)/L.Req
                        UL = -L.De*(1 + L.a1*y + L.a2*y**2 + L.a3*y**3 +
                                    L.a4*y**4 + L.a5*y**5) * np.exp(-L.a1*y) + L.Eref
                        return UL

                class complexa:
                    def UM_000(r):
                        UM_000 = (2*LC.UH(r)+LC.UL(r) + 2 *
                                  (LC.UTa(r)+LC.UTb(r)+LC.UX(r)))/9
                        return UM_000

                    def UM_202(r):
                        UM_202 = 2*(LC.UH(r)-LC.UL(r)+LC.UTa(r) -
                                    2*LC.UTb(r)+LC.UX(r))/(9*(5**(1/2)))
                        return UM_202

                    def UM_022(r):
                        UM_022 = 2*(LC.UH(r)-LC.UL(r)-2*LC.UTa(r) +
                                    LC.UTb(r)+LC.UX(r))/(9*(5**(1/2)))
                        return UM_022

                    def UM_220(r):
                        UM_220 = 2*(4*LC.UH(r)-LC.UL(r)-5*(LC.UTa(r) + 
                                    LC.UTb(r) + LC.UX(r))+12*LC.UZ(r))/(45*(5**(1/2)))
                        return UM_220

                    def UM_222(r):
                        UM_222 = ((2/7)**(1/2))*(13*LC.UH(r)-LC.UL(r)+7 *
                                                 (LC.UTa(r)+LC.UTb(r)-2*LC.UX(r))-12*LC.UZ(r))/45
                        return UM_222

                    def UM_224(r):
                        UM_224 = ((2/35)**(1/2)*8 * LC.UH(r)+LC.UL(r)+2*LC.UZ(r))/15
                        return UM_224

                def UM_FINAL(r, th_a, th_b, phi):  # potencial que vai ser usado na equação 3
                    UM_FINAL = (complexa.UM_000(r)+(5**(1/2))/4 * (3*(np.cos(2*th_a))+1) * complexa.UM_202(r) 
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
                    # print('um',UM_FINAL)
                    return UM_FINAL

            class drv:

                def d1_H(r):
                    d1_H = num_diff.derivative(
                        LC.UH, r, method='central', p=1e-8)
                    return d1_H

                def d1_L(r):
                    d1_L = num_diff.derivative(
                        LC.UL, r, method='central', p=1e-8)
                    return d1_L

                def d1_Ta(r):
                    d1_Ta = num_diff.derivative(
                        LC.UTa, r, method='central', p=1e-8)
                    return d1_Ta

                def d1_Tb(r):
                    d1_Tb = num_diff.derivative(
                        LC.UTb, r, method='central', p=1e-8)
                    return d1_Tb

                def d1_X(r):
                    d1_X = num_diff.derivative(
                        LC.UX, r, method='central', p=1e-8)
                    return d1_X

                def d1_Z(r):
                    d1_Z = num_diff.derivative(
                        LC.UZ, r, method='central', p=1e-8)
                    return d1_Z

                def d2_H(r):
                    d2_H = num_diff.derivative2(
                        LC.UH, r, method='central', p=1e-8)
                    return d2_H

                def d2_L(r):
                    d2_L = num_diff.derivative2(
                        LC.UL, r, method='central', p=1e-8)
                    return d2_L

                def d2_Ta(r):
                    d2_Ta = num_diff.derivative2(
                        LC.UTa, r, method='central', p=1e-8)
                    return d2_Ta

                def d2_Tb(r):
                    d2_Tb = num_diff.derivative2(
                        LC.UTb, r, method='central', p=1e-8)
                    return d2_Tb

                def d2_X(r):
                    d2_X = num_diff.derivative2(
                        LC.UX, r, method='central', p=1e-8)
                    return d2_X

                def d2_Z(r):
                    d2_Z = num_diff.derivative2(
                        LC.UZ, r, method='central', p=1e-8)
                    return d2_Z

                def d1r(r, th_a, th_b, phi):
                    d1r = (-1/18)*(1+3*np.cos(2*th_a))*(drv.d1_H(r) -
                                                        drv.d1_L(r)+drv.d1_Ta(r)-2*drv.d1_Tb(r)+drv.d1_X(r))
                    -(1/18)*(1+3*np.cos(2*th_b))*(drv.d1_H(r) -
                                                  drv.d1_L(r)-2*drv.d1_Ta(r)+drv.d1_Tb(r)+drv.d1_X(r))
                    +(1/9)*(2*drv.d1_H(r)+drv.d1_L(r)+2 *
                            (drv.d1_Ta(r)+drv.d1_Tb(r)+drv.d1_X(r)))
                    -(1/504)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))
                              - 3*(1-np.cos(2*th_a)) *
                              (1-np.cos(2*th_b))*np.cos(2*phi)
                              + 6*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(13*drv.d1_H(r)
                                                                              - drv.d1_L(r)+7*drv.d1_Ta(r)+7*drv.d1_Tb(r)-14*drv.d1_X(r)-12*drv.d1_Z(r))
                    +1/35*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+1/2*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)
                           - 8*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(drv.d1_H(r)+drv.d1_L(r)-2*drv.d1_Z(r))
                    +1/360*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)
                            + 12*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(4*drv.d1_H(r)-drv.d1_L(r)-5*drv.d1_Ta(r)-5*drv.d1_Tb(r)-5*drv.d1_X(r)+12*drv.d1_Z(r))
                    return d1r

                def d1th_a(r, th_a, th_b, phi):
                    d1th_a = 1/3*np.sin(2*th_a) * (LC.UH(r) -
                                                   LC.UL(r)+LC.UTa(r)-2*LC.UTb(r)+LC.UX(r))
                    -1/504*(-6*(1+3*np.cos(2*th_b)) * np.sin(2*th_a) - 6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
                            + 12*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
                    +1/35*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + (1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
                           - 16*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
                    +1/360*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + 6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)
                            + 24*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))

                    return d1th_a

                def d1th_b(r, th_a, th_b, phi):
                    d1th_b = 1/3*np.sin(2*th_b) * (LC.UH(r) -
                                                   LC.UL(r)-2*LC.UTa(r)+LC.UTb(r)+LC.UX(r))
                    -1/504*(12*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a) - 6*(1+3*np.cos(2*th_a)) * np.sin(2*th_b) - 6*(1-np.cos(
                        2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
                    +1/35*(-16*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a)
                           - 6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b) + (1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
                    +1/360*(24*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a) - 6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b)
                            + 6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b)) * (4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
                    return d1th_b

                def d1phi(r, th_a, th_b, phi):
                    d1phi = -1/504*(-6*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) + 6*(1-np.cos(2*th_a)) * (1-np.cos(
                        2*th_b)) * np.sin(2*phi)) * (13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
                    +1/35*(8*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - (1-np.cos(2*th_a))
                           * (1-np.cos(2*th_b)) * np.sin(2*phi)) * (LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
                    + 1/360*(-12*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - 6*(1-np.cos(2*th_a)) * (1-np.cos(
                        2*th_b)) * np.sin(2*phi))*(4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
                    return d1phi

                def d2r(r, th_a, th_b, phi):
                    d2r = -1/18*(1+3*np.cos(2*th_a)) * (drv.d2_H(r) -
                                                        drv.d2_L(r)+drv.d2_Ta(r)-2*drv.d2_Tb(r)+drv.d2_X(r))
                    -1/18*(1+3*np.cos(2*th_b)) * (drv.d2_H(r) -
                                                  drv.d2_L(r)-2*drv.d2_Ta(r)+drv.d2_Tb(r)+drv.d2_X(r))
                    +1/9 * (2*drv.d2_H(r)+drv.d2_L(r)+2 *
                            (drv.d2_Ta(r)+drv.d2_Tb(r)+drv.d2_X(r)))
                    -1/504*(1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b))
                    -3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi) + 6*np.cos(phi) * np.sin(2*th_a) * np.sin(
                        2*th_b) * (13*drv.d2_H(r)-drv.d2_L(r)+7*drv.d2_Ta(r)+7*drv.d2_Tb(r)-14*drv.d2_X(r)-12*drv.d2_Z(r))
                    +1/35*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 1/2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
                           - 8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (drv.d2_H(r)+drv.d2_L(r)-2*drv.d2_Z(r))
                    +1/360*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)
                            + 12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*drv.d2_H(r)-drv.d2_L(r)-5*drv.d2_Ta(r)-5*drv.d2_Tb(r)-5*drv.d2_X(r)+12*drv.d2_Z(r))
                    return d2r

                def d2th_a(r, th_a, th_b, phi):
                    d2th_a = 2/3*np.cos(2*th_a) * (LC.UH(r) -
                                                   LC.UL(r)+LC.UTa(r)-2*LC.UTb(r)+LC.UX(r))
                    + 1/504*(12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
                             + 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*LC.UH(r)-LC.UL(r)+7*LC.UTa(r)+7*LC.UTb(r)-14*LC.UX(r)-12*LC.UZ(r))
                    + 1/35*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 2*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
                            + 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (LC.UH(r)+LC.UL(r)-2*LC.UZ(r))
                    + 1/360*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)
                             - 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*LC.UH(r)-LC.UL(r)-5*LC.UTa(r)-5*LC.UTb(r)-5*LC.UX(r)+12*LC.UZ(r))
                    return d2th_a

                def d2th_b(r, th_a, th_b, phi):
                    d2th_b = 2/3*np.cos(2*th_b) * (LC.UH(r) -
                                                   LC.UL(r)-2*LC.UTa(r)+LC.UTb(r)+LC.UX(r))
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
                def BH2_ref(T):  # from 60K to 500K
                    return 1.7472*10 - (1.2926 * 10**2)/(T) - (2.6988 * 10**5)/(T)**2 + (8.0282 * 10**6)/(T)**3

                def BCl2_ref(T):  # from  300K  to 1070?
                    return 1.3171*10 + 5.1055*10**3/T - 2.9404*10**7/T**2 - 5.0383*10**8/T**3

                def BF2_ref(T):  # from  85K  to 295?
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

            Reqs = np.array([H.Req, X.Req, Z.Req, Ta.Req, Tb.Req, L.Req])

            media_Reqs = mean(Reqs)
            if gas == 'F2F2':
                #st.write('f2f2')
                m1 = massaF
                m2 = massaF
                m3 = massaF
                m4 = massaF

                Req1 = 1.415962
                Req2 = 1.415962  # Cl  2.007880 # F = 1.415962 #usar dados CBS da tabela

                T_ref = [80, 90, 100, 110, 120, 160, 200, 250, 300]  # errado
                B_ref = [-50.13, -50.13, -33.91, -27.86,  -22.83,  -
                         9.18,  -1.30,  4.70, 8.49]  # errado, sao do h2f2
                
                mi = (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
                #mi1 = (m1 + m2) / 4 * Req1**2 / 6.02252
                #mi2 = (m3 + m4) / 4 * Req2**2 / 6.02252
                
                
                
                r1 = (m2 * media_Reqs) / (m1 - m2)
                r2 = media_Reqs - r1
                r3 = (m4 * media_Reqs) / (m3 - m4)
                r4 = media_Reqs - r3



                I1 = m1 * r1**2 + m2 * r2**2
                I2 =  m3 * r3**2 + m4 * r4**2


                mi1 = (m1 * m2)/(m1 + m2)
                mi2 = (m3 * m4)/(m3 + m4)


                Ma = m1 + m2
                Mb = m3 + m4
                w1 = 0.053
                w2 = 0.053
                rho_1 = 0.02715
                rho_2 = 0.02715
                zc1 = 0.287
                zc2 = 0.287
                x1 = 0.5
                x2 = 0.5

                def B_A2_ref(T):
                    return B_references.BF2_ref(T)

                def B_B2_ref(T):
                    return B_references.BF2_ref(T)

            elif gas == 'H2F2':
                #st.write('h2f2')
                m1 = massaH
                m2 = massaH
                m3 = massaF
                m4 = massaF

                Req1 = 0.744013
                Req2 = 1.415962  # Cl  2.007880 # F = 1.415962 #usar dados CBS da tabela
                T_ref = [80, 90, 100, 110, 120, 160, 200, 250, 300]
                B_ref = [-50.13, -50.13, -33.91, -27.86,  -
                         22.83,  -9.18,  -1.30,  4.70, 8.49]

                mi = (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
                #mi1 = (m1 + m2) / 4 * Req1**2 / 6.02252
                #mi2 = (m3 + m4) / 4 * Req2**2 / 6.02252

                r1 = (m2 * media_Reqs) / (m1 - m2)
                r2 = media_Reqs - r1
                r3 = (m4 * media_Reqs) / (m3 - m4)
                r4 = media_Reqs - r3



                I1 = m1 * r1**2 + m2 * r2**2
                I2 =  m3 * r3**2 + m4 * r4**2


                mi1 = (m1 * m2)/(m1 + m2)
                mi2 = (m3 * m4)/(m3 + m4)                


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

            elif gas == 'H2H2':
                #st.write('h2h2')

                m1 = massaH
                m2 = massaH
                m3 = massaH
                m4 = massaH
                Req1 = 0.744013
                Req2 = 0.744013

                mi = (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
                #mi1 = (m1 + m2) / 4 * Req1**2 / 6.02252
                #mi2 = (m3 + m4) / 4 * Req2**2 / 6.02252

                r1 = (m2 * media_Reqs) / (m1 - m2)
                r2 = media_Reqs - r1
                r3 = (m4 * media_Reqs) / (m3 - m4)
                r4 = media_Reqs - r3



                I1 = m1 * r1**2 + m2 * r2**2
                I2 =  m3 * r3**2 + m4 * r4**2


                mi1 = (m1 * m2)/(m1 + m2)
                mi2 = (m3 * m4)/(m3 + m4)

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

            elif gas == 'H2Br2':
                #st.write('h2Br2')
                m1 = massaH
                m2 = massaH
                m3 = massaBr
                m4 = massaBr
                B_ref = []
                T_ref = []
                Req1 = 0
                Req2 = 0

                mi = (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
                #mi1 = (m1 + m2) / 4 * Req1**2 / 6.02252
                #mi2 = (m3 + m4) / 4 * Req2**2 / 6.02252

                r1 = (m2 * media_Reqs) / (m1 - m2)
                r2 = media_Reqs - r1
                r3 = (m4 * media_Reqs) / (m3 - m4)
                r4 = media_Reqs - r3



                I1 = m1 * r1**2 + m2 * r2**2
                I2 =  m3 * r3**2 + m4 * r4**2


                mi1 = (m1 * m2)/(m1 + m2)
                mi2 = (m3 * m4)/(m3 + m4)

                Ma = m1 + m2  # verify
                Mb = m3 + m4  # esses dados sao do h2br2
                w1 = -0.216
                w2 = 0.286
                rho_1 = 0.026667
                rho_2 = 0.05538
                zc1 = 0.303
                zc2 = 0.126
                x1 = 0.5
                x2 = 0.5

            else:
                #st.write('h2cl2')
                m1 = massaH
                m2 = massaH
                m3 = massaCl
                m4 = massaCl
                B_ref = [-66.0103, -52.9951, -42.88, -34.849, -28.3552, -23.0219,  -18.5833,  -14.8471, -11.6713, -8.94904,  -6.59832, -4.55533, -2.76969, -1.20109,  0.18304, 1.40931,
                         2.49968,  3.47237,  4.34266, 5.12343, 5.82562, 6.45854,  7.0302, 7.54749, 8.01638,  8.44207, 8.82907, 9.18137, 9.50244,  9.79537,  10.0629,  10.3074,  10.531]
                T_ref = [300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700,
                         725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000, 1025,  1050,  1075, 1100, ]
                Req1 = 0.744013
                Req2 = 2.007880

                mi = (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
                #mi1 = (m1 + m2) / 2 * Req1**2 #/ 6.02252
                #mi2 = (m3 + m4) / 2 * Req2**2 #/ 6.02252

                r1 = (m2 * media_Reqs) / (m1 - m2)
                r2 = media_Reqs - r1
                r3 = (m4 * media_Reqs) / (m3 - m4)
                r4 = media_Reqs - r3



                I1 = m1 * r1**2 + m2 * r2**2
                I2 =  m3 * r3**2 + m4 * r4**2


                mi1 = (m1 * m2)/(m1 + m2)
                mi2 = (m3 * m4)/(m3 + m4)

                Ma = m1 + m2  # verify
                Mb = m3 + m4  # esses dados sao do h2cl2
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
                
            #st.write('media_Reqs = ', media_Reqs, 'r1 = ',r1, 'r2 = ',r2, 'r3 = ', r3, 'r4 = ',r4 ,  'I1 = ', I1,  'I2 = ', I2,  'mi1 = ', mi1,  'mi2 = ', mi2)

            B_virial_state_ref = []
            T_state = []

            Reqs = np.array([H.Req, X.Req, Z.Req, Ta.Req, Tb.Req, L.Req])

            lims = mean(Reqs)
            lim_inf = lims/4
            lim_sup = 10*lims/2
            r0 = 4.22  # mean(Reqs)

            for T in range(60, 500, 5):

                B_cruzado = (0.25*(22.1*Mb + Ma*(22.1 + Mb*T)) * (0.083 + 0.0695*w1 + 0.0695*w2) * (
                    rho_1**(1/3) + rho_2**(1/3))**3) / ((10.9*Ma + 10.9*Mb + Ma*Mb*T)*(zc1 + zc2) * rho_1 * rho_2)

                B_virial_state = 2*x1*x2*B_cruzado + \
                    x2**2*B_A2_ref(T) + x1**2*B_B2_ref(T)

                B_virial_state_ref.append(B_virial_state)
                T_state.append(T)


            def integrand_vegas(x):
                r = x[0]
                th_a = x[1]
                th_b = x[2]
                phi = x[3]
                F = UM_FINAL(r, th_a, th_b, phi)*constants.rk/T

                # print('f',F,x)

                # *10e-24 #24 eh o certo de ang^3 pra cm ^3
                return constants.N_A/4 * np.sin(th_a)  * np.sin(th_b)*(r**2)*(1-np.exp(-F))

            def integrand_c1(x):  # primeira correção quantica
                r = x[0]
                th_a = x[1]
                th_b = x[2]
                phi = x[3]
                F = UM_FINAL(r, th_a, th_b, phi)*constants.rk/T

                d_1r = drv.d1r(r, th_a, th_b, phi)
                #c1 = N_A * h**2 * rk**3 / (48 * mi * T**3) * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * d_1r**2
                c1 =  -(constants.N_A * constants.h**2 * constants.rk**3)/(96*mi*T**3) * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * d_1r**2 * r**2

                return c1

            def integrand_c2(x):  # segunda correção quantica
                r = x[0]
                th_a = x[1]
                th_b = x[2]
                phi = x[3]
                F = UM_FINAL(r, th_a, th_b, phi)*constants.rk/T

                d_1th_a = drv.d1th_a(r, th_a, th_b, phi)
                d_1th_b = drv.d1th_b(r, th_a, th_b, phi)
                d_2th_a = drv.d2th_a(r, th_a, th_b, phi)
                d_2th_b = drv.d2th_b(r, th_a, th_b, phi)
                d_2phi = drv.d2phi(r, th_a, th_b, phi)

                VL1 = (constants.h**2)/(2*I1) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/np.sin(th_a)**2 * d_2phi)
                VL2 = (constants.h**2)/(2*I2) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/np.sin(th_b)**2 * d_2phi)

                
                #d_1th_a * np.cos(th_a)/np.sin(th_a) - \
                #    d_2th_a - d_2phi/(np.sin(th_a))**2
                #VL2 = (constants.h**2)/(2*mi1*r**2) * d_1th_b * np.cos(th_b)/np.sin(th_b) - \
                #    d_2th_b - d_2phi/(np.sin(th_b))**2

                c2 = -(constants.N_A * constants.rk**3)/(48*T**3) * np.exp(-F) * (VL1 + VL2) * r**2

                # c2a =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/(np.power(np.sin(th_a),2)) * d_2phi) * r**2  #talvez sem esse rk/T
                # -np.pi/12 * N_A * h**2/mi1 *
                # c2b =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/(np.power(np.sin(th_b),2)) * d_2phi) * r**2  #talvez sem esse rk/T

                #c2 = -np.pi/12 * N_A * h**2 * (c2a/mi1 +  c2b/mi2) * (rk/T)**2
                return c2

            def integrand_c3(x):  # segunda correção quantica
                r = x[0]
                th_a = x[1]
                th_b = x[2]
                phi = x[3]
                F = UM_FINAL(r, th_a, th_b, phi)*constants.rk/T

                d_1th_a = drv.d1th_a(r, th_a, th_b, phi)
                d_1th_b = drv.d1th_b(r, th_a, th_b, phi)
                d_2th_a = drv.d2th_a(r, th_a, th_b, phi)
                d_2th_b = drv.d2th_b(r, th_a, th_b, phi)
                d_2phi = drv.d2phi(r, th_a, th_b, phi)

                VL1 = (constants.h**2)/(2*mi1*r**2) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/np.sin(th_a)**2 * d_2phi)
                VL2 = (constants.h**2)/(2*mi2*r**2) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/np.sin(th_b)**2 * d_2phi)

                
                #d_1th_a * np.cos(th_a)/np.sin(th_a) - \
                #    d_2th_a - d_2phi/(np.sin(th_a))**2
                #VL2 = (constants.h**2)/(2*mi1*r**2) * d_1th_b * np.cos(th_b)/np.sin(th_b) - \
                #    d_2th_b - d_2phi/(np.sin(th_b))**2

                c3 = -(constants.N_A * constants.rk**3)/(48*T**3) * np.exp(-F) * (VL1 + VL2) * r**2

                # c2a =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/(np.power(np.sin(th_a),2)) * d_2phi) * r**2  #talvez sem esse rk/T
                # -np.pi/12 * N_A * h**2/mi1 *
                # c2b =  np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/(np.power(np.sin(th_b),2)) * d_2phi) * r**2  #talvez sem esse rk/T

                #c2 = -np.pi/12 * N_A * h**2 * (c2a/mi1 +  c2b/mi2) * (rk/T)**2
                return c3

            def integrand_c4(x):  # quarta correção quantica
                r = x[0]
                th_a = x[1]
                th_b = x[2]
                phi = x[3]
                F = UM_FINAL(r, th_a, th_b, phi)*constants.rk/T

                d_1r = drv.d1r(r, th_a, th_b, phi)
                d_2r = drv.d2r(r, th_a, th_b, phi)
                d_1th_a = drv.d1th_a(r, th_a, th_b, phi)
                d_1th_b = drv.d1th_b(r, th_a, th_b, phi)
                d_2th_a = drv.d2th_a(r, th_a, th_b, phi)
                d_2th_b = drv.d2th_b(r, th_a, th_b, phi)
                d_2phi = drv.d2phi(r, th_a, th_b, phi)



                c4 = -(constants.N_A * constants.h**4 * constants.rk**4)/(1920 * mi**2 * T**4) * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * \
                (d_2r**2  + 2/r**2 * d_1r**2 + (10 * constants.rk)/(9*T*r) * d_1r**3 - (5* constants.rk**2)/(36*T**2) * d_1r**4) * r**2                                      



                #c4a = -N_A * h**4 * rk**4 / (1920 * mi**2 * T**4) * np.sin(th_a) * np.sin(th_b) * np.exp(-F)
                #c4b = (d_2r**2 + 2/r**2*(d_1r**2) + 10 * rk/(9 * T * r) * d_1r**3 - 5 * rk**2/(36 * T**2) * d_1r**4) * r**2
                #c4aa = c4a * c4b
                # c4p = np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (d_2r**2 + 2/r**2*(d_1r**2) + 10 * rk/(9 * T * r) * d_1r**3 - 5 * rk**2/(36 * T**2) * d_1r**4) * r**2  #talvez sem rk/T
                #c4 = -np.pi/12 * N_A * h**2/mi * c4p *(rk/T)**2
                return c4

            integ = vegas.Integrator(
                [[0, 10], [0, np.pi], [0, np.pi], [0, 2*np.pi]])

            B_clas = []
            Tstr = []

            B_plus_c1 = []
            B_main = []
            B_c1 = []
            B_c2 = []
            B_c3 = []
            B_c4 = []
            B_correcoes = []
            B_plus_all_except_c2 = []

            T = step
            integ(integrand_vegas, nitn=10, neval=1000)

            integ(integrand_c1, nitn=10, neval=1000)

            integ(integrand_c2, nitn=10, neval=1000)
            integ(integrand_c3, nitn=10, neval=1000)
            integ(integrand_c4, nitn=10, neval=1000)

            result = integ(integrand_vegas, nitn=10, neval=1000)

            # st.write(result.summary())

            st.write('result of classic virial = ', result, 'for T = ', T)
            
            B_clas.append(result.mean)
            Tstr.append(T)

            result_c1 = integ(integrand_c1, nitn=10, neval=1000)

            result_c2 = integ(integrand_c2, nitn=10, neval=1000)
            result_c3 = integ(integrand_c3, nitn=10, neval=1000)
            result_c4 = integ(integrand_c4, nitn=10, neval=1000)

            #st.write('result VEGAS c1 = ', result_c1)
            #st.write('result VEGAS c2 = ', result_c2)
            #st.write('result VEGAS c3 = ', result_c3)
            #st.write('result VEGAS c4 =', result_c4)

            st.write('result of final virial =', result.mean + result_c1.mean +
                        result_c2.mean + result_c3.mean + result_c4.mean, 'for T = ', T)

            B_main.append(result.mean + result_c1.mean +
                            result_c2.mean + result_c3.mean + result_c4.mean)

            B_c1.append(result_c1.mean)
            B_c2.append(result_c2.mean)
            B_c3.append(result_c3.mean)
            B_c4.append(result_c4.mean)
            B_correcoes.append(
                +result_c2.mean + result_c3.mean + result_c1.mean + result_c4.mean)

            B_plus_all_except_c2.append(
                +result.mean + result_c1.mean + result_c3.mean + result_c4.mean)
            


            data_virial = data_virial.append({'Temperature':T,'Classical Virial Coefficient':result.mean, 'First Virial Correction':result_c1.mean, 
                'Second Virial Correction':result_c2.mean, 'Third Virial Correction':result_c3.mean, 'Fourth Virial Correction':result_c4.mean,},ignore_index=True)
            

            st.write(data_virial)
            r = np.linspace(0, 10, 100)
            th_a = np.linspace(0, np.pi, 100)
            th_b = np.linspace(0, np.pi, 100)
            phi = np.linspace(0, 2*np.pi, 100)

            #fig = make_subplots(rows=1, cols=3, subplot_titles=("(a)energias simples", "(b)energias complexas", "(c) Segundo Coeficiente Virial(B) em função da Temperatura(T)"))
            fig1 = go.Figure()

            fig1.add_trace(
                go.Scatter(x=r, y=LC.UH(r), name='UH')
            )
            fig1.add_trace(
                go.Scatter(x=r, y=LC.UX(r), name='UX')
            )
            fig1.add_trace(
                go.Scatter(x=r, y=LC.UZ(r), name='UZ')
            )
            fig1.add_trace(
                go.Scatter(x=r, y=LC.UTa(r), name='UTa')
            )
            fig1.add_trace(
                go.Scatter(x=r, y=LC.UTb(r), name='UTb')
            )
            fig1.add_trace(
                go.Scatter(x=r, y=LC.UL(r), name='UL')
            )
            # ************************ax1.set_ylabel(r'$U(r)[eV]$')
            # ************************ax1.set_xlabel(r'$r[Å]$', labelpad=1)

            # ax1.set_xlim([2, 8]) #1.2 B
            # ax1.set_ylim([-40, 100]) #inferior : B(T=10) * 1.1 superior: B(T=500)* 1.2

            # ax2.set_xlim([2, 8]) #1.2 B
            #ax2.set_ylim([-40, 100])
            # ax3.set_xlim([2, 10]) #1.2 B
            #ax3.set_ylim([-1000, 1006])
            # ax4.set_xlim([2, 10]) #1.2 B
            #ax4.set_ylim([-1000, 1006])

            fig2 = go.Figure()
            fig2.add_trace(
                go.Scatter(x=r, y=complexa.UM_000(r), name='UM_000')
            )
            fig2.add_trace(
                go.Scatter(x=r, y=complexa.UM_202(r), name='UM_202')
            )
            fig2.add_trace(
                go.Scatter(x=r, y=complexa.UM_220(r), name='UM_220')
            )
            fig2.add_trace(
                go.Scatter(x=r, y=complexa.UM_022(r), name='UM_022')
            )
            fig2.add_trace(
                go.Scatter(x=r, y=complexa.UM_222(r), name='UM_222')
            )
            fig2.add_trace(
                go.Scatter(x=r, y=complexa.UM_224(r), name='UM_224')
            )

            # ax2.set_xlim([2, 9]) #1.2 B
            # ax2.set_ylim([-0.4, 2]) #inferior : B(T=10) * 1.1 superior: B(T=500)* 1.2
            # ******************ax2.set_ylabel(r'$U(r)[eV]$')
            # *****************ax2.set_xlabel(r'$r[Å]$', labelpad=1)



            # ax3.plot(T_ref, B_ref,   color='b',label = 'ref') #Breferencia

            #ax3.plot(Tstr, B_main,   color='yellow',label = 'Bmain')

            # ****************** ax3.set_ylabel(r'$B[cm^3/mol]$')
            # ***************ax3.set_xlabel(r'$T[K]$', labelpad=1)
        
            #plt.subplot(Tstr, Bstr)
            #plt.scatter(Tstr, Bstr)
            #plt.title(f'Gráficos (a) do , name. You are age.')

            fig1.update_xaxes(range=[0, 10])
            fig1.update_yaxes(range=[-3*H.De, 3*H.De])

            fig2.update_xaxes(range=[0, 10])
            fig2.update_yaxes(range=[-3*H.De, 5*H.De])

            #fig3.update_xaxes(range=[0, 1000])
            #fig3.update_yaxes(range=[-1000, 1000000000])

            #fig.update_yaxes(title_text="yaxis 3 title", showgrid=False)

            #fig.update_xaxes( range=[0, 5], showgrid=False)

            fig1.update_layout(height=600, width=800,
                               title_text="Graph of Leading Configurations per Distance", xaxis_title=r"$r[A]$",
                               yaxis_title='U[eV]'
                                )
            st.plotly_chart(fig1, use_container_width=True)

            fig2.update_layout(height=600, width=800,
                               title_text="Graph of Complex Energies per Distance", xaxis_title=r"$r[A]$",
                               yaxis_title='U[eV]')
            st.plotly_chart(fig2, use_container_width=True)





            # st.pyplot(f)
            #st.pyplot(f, ax1, ax2,ax3)
            st.write('evaluation time {}'.format(time.strftime(
                "%H:%M:%S", time.gmtime(time.time()-start))))
            st.write('number of temperatures:', len(Tstr))

            st.success('Calculations finished!', icon="✅")
        elif gas_type == 'AB':
            st.write('If you wish to calculate the second virial coefficient of a diatomic molecule, please reffer to the page "diatomic molecules AB')





    else:
        st.info('☝️ Upload a .dat file')
