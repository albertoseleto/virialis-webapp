# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:00:27 2023

@author: Alberto
"""
#########################
#LIBRARY IMPORTS
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
import logging


from numeric import num_diff
from constant import constants
from collections import namedtuple

###############################

# TIME MEASURE

# DATA ENTRY



logging.basicConfig(filename='log_information.log', encoding='utf-8', level=logging.DEBUG)


# Define a named tuple to represent the return values from function_A

st.title('Second Virial Coefficient Calculator for A2B2 Molecule')

uploaded_file = st.file_uploader("upload a file")


potential = st.selectbox(
    'What potential energy do you want to use?',
    ('Rydberg Potential', 'Improved Leonard-Jonnes Potential', 'both'))


step = st.selectbox(
    'What is the Temperature step that you want to use?',
    (100, 50, 25, 200, 300))

gas_type = st.selectbox(
    'What gas structure are you studying?',
    ('A2B2', 'AB',))
if potential == 'both':
    uploaded_file_2 = st.file_uploader("upload the rydberg file")

if gas_type == 'A2B2':

    monomero1 = st.selectbox(
        'What is your first monomer?',
        ('H2', 'F2', 'Cl2', 'Br2'))
    monomero2 = st.selectbox(
        'What is your second monomer?',
        ('H2', 'F2', 'Cl2', 'Br2'))
    gas = monomero1 + monomero2


st.write('The calculus will be made using ', potential,
          'using a temperature step of', step)

ILJ_result = namedtuple('ILJresult', ["Reqs", "LC_lst", "harm_esf_lst", "Tstr", "B_clas", "De_lim"])

def ILJ():
    st.write('LEONARD-JONNES IS IN CONSTRUCTION...')

    data = pd.read_csv(uploaded_file, sep="\s+", header=None)
    data.columns = ["alpha", "beta", "mp", "De", "Req" ]
    st.subheader('DataFrame')
    st.write(data)


    #dataframes de input
    h_alpha = data.loc[0, 'alpha']
    h_beta = data.loc[0, 'beta']
    h_mp = data.loc[0, 'mp']
    h_De = data.loc[0, 'De'] 
    h_Req = data.loc[0, 'Req']
        

    x_alpha = data.loc[1, 'alpha']
    x_beta = data.loc[1, 'beta']
    x_mp = data.loc[1, 'mp']
    x_De = data.loc[1, 'De'] 
    x_Req = data.loc[1, 'Req']
        

    z_alpha = data.loc[2, 'alpha']
    z_beta = data.loc[2, 'beta']
    z_mp = data.loc[2, 'mp']
    z_De = data.loc[2, 'De'] 
    z_Req = data.loc[2, 'Req']
        

    ta_alpha = data.loc[3, 'alpha']
    ta_beta = data.loc[3, 'beta']
    ta_mp = data.loc[3, 'mp']
    ta_De = data.loc[3, 'De'] 
    ta_Req = data.loc[3, 'Req']
    

    tb_alpha = data.loc[4, 'alpha']
    tb_beta = data.loc[4, 'beta']
    tb_mp = data.loc[4, 'mp']
    tb_De = data.loc[4, 'De'] 
    tb_Req = data.loc[4, 'Req']
        

    l_alpha = data.loc[5, 'alpha']
    l_beta = data.loc[5, 'beta']
    l_mp = data.loc[5, 'mp']
    l_De = data.loc[5, 'De'] 
    l_Req = data.loc[5, 'Req']
        


    #def de funcoes
    class LC1: 
        def UH(r):
            n = h_beta + h_alpha * (r / h_Req) ** 2

            return  h_De * ((h_mp/(n - h_mp) * (h_Req/r) ** n) - (n/(n - h_mp) * (h_Req/r) ** h_mp))

        def UX(r):

            n = x_beta + x_alpha * (r / x_Req) ** 2
            return x_De * ((x_mp/(n - x_mp) * (x_Req/r) ** n) - (n/(n - x_mp) * (x_Req/r) ** x_mp))

        def UZ(r):

            n = z_beta + z_alpha * (r / z_Req) ** 2
            return z_De * ((z_mp/(n - z_mp) * (z_Req/r) ** n) - (n/(n - z_mp) * (z_Req/r) ** z_mp))

        def UTa(r):
            n = ta_beta + ta_alpha * (r / ta_Req) ** 2
            return ta_De * ((ta_mp/(n - ta_mp) * (ta_Req/r) ** n) - (n/(n - ta_mp) * (ta_Req/r) ** ta_mp))

        def UTb(r):
            n = tb_beta + tb_alpha * (r / tb_Req) ** 2
            return tb_De * ((tb_mp/(n - tb_mp) * (tb_Req/r) ** n) - (n/(n - tb_mp) * (tb_Req/r) ** tb_mp))

        def UL(r):
            n = l_beta + l_alpha * (r / l_Req) ** 2
            return l_De * ((l_mp/(n - l_mp) * (l_Req/r) ** n) - (n/(n - l_mp) * (l_Req/r) ** l_mp))
        '''
        def USa(r):
            n = Sa.beta + Sa.alpha * (r / Sa.Req) ** 2
            return Sa.De * ((Sa.mp/(n - Sa.mp) * (Sa.Req/r) ** n) - (n/(n - Sa.mp) * (Sa.Req/r) ** Sa.mp))

        def USb(r):
            n = Sb.beta + Sb.alpha * (r / Sb.Req) ** 2
            return Sb.De * ((Sb.mp/(n - Sb.mp) * (Sb.Req/r) ** n) - (n/(n - Sb.mp) * (Sb.Req/r) ** Sb.mp))
        '''

    class complexa1:
        def UM_000(r):
            UM_000 = (2*LC1.UH(r)+LC1.UL(r) + 2 *
                        (LC1.UTa(r)+LC1.UTb(r)+LC1.UX(r)))/9
            return UM_000

        def UM_202(r):
            UM_202 = 2*(LC1.UH(r)-LC1.UL(r)+LC1.UTa(r) -
                        2*LC1.UTb(r)+LC1.UX(r))/(9*(5**(1/2)))
            return UM_202

        def UM_022(r):
            UM_022 = 2*(LC1.UH(r)-LC1.UL(r)-2*LC1.UTa(r) +
                        LC1.UTb(r)+LC1.UX(r))/(9*(5**(1/2)))
            return UM_022

        def UM_220(r):
            UM_220 = 2*(4*LC1.UH(r)-LC1.UL(r)-5*(LC1.UTa(r) + 
                        LC1.UTb(r) + LC1.UX(r))+12*LC1.UZ(r))/(45*(5**(1/2)))
            return UM_220

        def UM_222(r):
            UM_222 = ((2/7)**(1/2))*(13*LC1.UH(r)-LC1.UL(r)+7 *
                                        (LC1.UTa(r)+LC1.UTb(r)-2*LC1.UX(r))-12*LC1.UZ(r))/45
            return UM_222

        def UM_224(r):
            UM_224 = ((2/35)**(1/2)*8 * LC1.UH(r)+LC1.UL(r)+2*LC1.UZ(r))/15
            
            return UM_224

    def UM_FINAL1(r, th_a, th_b, phi):  # potencial que vai ser usado na equação 3
        UM_FINAL = (complexa1.UM_000(r)+(5**(1/2))/4 * (3*(np.cos(2*th_a))+1) * complexa1.UM_202(r) 
            + (5**(1/2))/4 * (3*(np.cos(2*th_b))+1) * complexa1.UM_022(r) 
            + (5**(1/2))/16 * complexa1.UM_220(r) * ((3*(np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
            + 12*np.sin(2*th_a)* np.sin(2*th_b)*np.cos(phi)
            + 3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)) 
            - 14**(1/2)*5/112 * complexa1.UM_222(r) * ((3*np.cos(2*th_a))+1) * (3*(np.cos(2*th_b))+1) 
            + 6*np.sin(2*th_a)*np.sin(2*th_b) * np.cos(phi) 
            - 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b))*np.cos(2*phi)
            + (3*(70)**(1/2))/112 * complexa1.UM_224(r) * ((3*(np.cos(2*th_a))+1)*(3*(np.cos(2*th_b))+1) 
            - 8* np.sin(2*th_a) * np.sin(2*th_b)*np.cos(phi)+((1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi))/2))
        #print('r:',r ,'000:',UM_000(r),'202:',UM_202(r),'022:',UM_022(r),'220:',UM_220(r),'222:',UM_222(r),'224:',UM_224(r))
        #st.write('um',UM_FINAL)
        

        return UM_FINAL
    

    #funcao final a ser integrada
    last_F_container = [-1e100]
    def integrand_vegas(x):
        last_F = last_F_container[0]
        r = x[0]
        th_a = x[1]
        th_b = x[2]
        phi = x[3]
        F = UM_FINAL1(r, th_a, th_b, phi)*constants.rk/T
        if F < -2:
            F = -2


        
        ff = constants.N_A/4 * np.sin(th_a)  * np.sin(th_b)*(r**2)*(1-np.exp(-F))
      
        #if F < -4000:
        #    logging.debug(f'ff value: {ff}, F value: {F}, r: {r}, th_a: {th_a}, th_b: {th_b}, phi: {phi}')
        #    return -1e30
        #else:
        #    return ff
        

        if ff < -1e4:
            if F >= last_F:
                logging.debug(f'ff value: {ff}, F value: {F}, r: {r}, th_a: {th_a}, th_b: {th_b}, phi: {phi}')
                last_F_container[0] = F
            return ff
        else:
            return ff


    Reqs = h_Req, x_Req, z_Req, ta_Req, tb_Req, l_Req
    media_Reqs = mean(Reqs)
    st.write(media_Reqs)

    #integrador e def de limites de integracao
    integ = vegas.Integrator(
        [[0.2*media_Reqs, 2*media_Reqs], [0, np.pi], [0, np.pi], [0, 2*np.pi]])
    #R0_inf=2.328d0  possivelmente nao ideal, manter 0  
    #RN_sup=23.813d0

    B_clas = []
    Tstr = []

    B_main = []
    
    for T in range(50, 1000, step):
        

        
        result = integ(integrand_vegas, nitn=10, neval=10000)

        # st.write(result.summary())

        st.write('result of classic virial = ', result, 'for T = ', T)
        
        B_clas.append(result.mean)
        Tstr.append(T)



    st.write('result of final virial =', result.mean, 'for T = ', T)

    B_main.append(result.mean )



    r = np.linspace(0, 10, 100)
    th_a = np.linspace(0, np.pi, 100)
    th_b = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)

    LC_lst = [LC1.UH(r), LC1.UX(r),LC1.UZ(r), LC1.UTa(r), LC1.UTb(r),LC1.UL(r)]

    harm_esf_lst = [complexa1.UM_000(r),complexa1.UM_202(r), complexa1.UM_022(r), complexa1.UM_220(r),complexa1.UM_222(r), complexa1.UM_224(r)]
    De_lim = h_De
    return ILJ_result(Reqs, LC_lst, harm_esf_lst, Tstr, B_clas, De_lim)


RYD_result = namedtuple('RYDresult', ["Reqs", "LC_lst", "harm_esf_lst", "Tstr", "B_clas", "De_lim"])

def ryd():
    if potential == "both":
        
        data = pd.read_csv(uploaded_file_2, sep="\s+", header=None)
    else: 
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
    
    #funcao final a ser integrada
    last_F_container = [0]
    def integrand_vegas(x):
        last_F = last_F_container[0]
        r = x[0]
        th_a = x[1]
        th_b = x[2]
        phi = x[3]
        F = UM_FINAL(r, th_a, th_b, phi)*constants.rk/T

        if F < -2:
            F = -2

        ff = constants.N_A/4 * np.sin(th_a)  * np.sin(th_b)*(r**2)*(1-np.exp(-F))
      
        #if F < -4000:
        #    logging.debug(f'ff value: {ff}, F value: {F}, r: {r}, th_a: {th_a}, th_b: {th_b}, phi: {phi}')
        #    return -1e30
        #else:
        #    return ff
        

    
        if F <= last_F:
            logging.debug(f'ff value: {ff}, F value: {F}, r: {r}, th_a: {th_a}, th_b: {th_b}, phi: {phi}')
            last_F_container[0] = F

 
        return ff
    
    Reqs = np.array([H.Req, X.Req, Z.Req, Ta.Req, Ta.Req, L.Req])
    media_Reqs = mean(Reqs)
    st.write(media_Reqs)

    #integrador e def de limites de integracao
    integ = vegas.Integrator(
        [[0.2*media_Reqs, 2*media_Reqs], [0, np.pi], [0, np.pi], [0, 2*np.pi]])

    B_clas = []
    Tstr = []

    B_main = []
    

    for T in range(50, 1000, step):

        integ(integrand_vegas, nitn=10, neval=10000)

        
        result = integ(integrand_vegas, nitn=10, neval=10000)

        # st.write(result.summary())

        st.write('result of classic virial = ', result, 'for T = ', T)
        
        B_clas.append(result.mean)
        Tstr.append(T)



    st.write('result of final virial =', result.mean, 'for T = ', T)

    B_main.append(result.mean )



    r = np.linspace(0, 10, 100)
    th_a = np.linspace(0, np.pi, 100)
    th_b = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)



    LC_lst = [LC.UH(r), LC.UX(r),LC.UZ(r), LC.UTa(r), LC.UTb(r),LC.UL(r)]

    harm_esf_lst = [complexa.UM_000(r),complexa.UM_202(r), complexa.UM_022(r), complexa.UM_220(r),complexa.UM_222(r), complexa.UM_224(r)]
    De_lim = H.De
    return RYD_result(Reqs, LC_lst, harm_esf_lst, Tstr, B_clas, De_lim)
def ref_making(Reqs):

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
    

    media_Reqs = mean(Reqs)
    massaCl = 34.96885
    if gas == 'H2F2':
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
        
    elif gas == 'F2F2':
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
                    725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000, 1025,  1050,  1075, 1100 ]
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
    return B_virial_state_ref, T_state

            
def graph_gen(LC_lst, harm_esf_lst, B_virial_state_ref, T_state, B_clas, Tstr, De_lim, 
              LC_lst_2=None, harm_esf_lst_2=None, B_clas_2=None, Tstr_2=None, De_lim_2=None):

    r = np.linspace(0, 10, 100)
    th_a = np.linspace(0, np.pi, 100)
    th_b = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)

    fig1 = go.Figure()

    fig1.add_trace(
        go.Scatter(x=r, y=LC_lst[0], name='UH')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC_lst[1], name='UX')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC_lst[2], name='UZ')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC_lst[3], name='UTa')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC_lst[4], name='UTb')
    )
    fig1.add_trace(
        go.Scatter(x=r, y=LC_lst[5], name='UL')
    )


    fig2 = go.Figure()
    fig2.add_trace(
        go.Scatter(x=r, y=harm_esf_lst[0], name='UM_000')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=harm_esf_lst[1], name='UM_202')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=harm_esf_lst[2], name='UM_220')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=harm_esf_lst[3], name='UM_022')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=harm_esf_lst[4], name='UM_222')
    )
    fig2.add_trace(
        go.Scatter(x=r, y=harm_esf_lst[5], name='UM_224')
    )


    fig3 = go.Figure()



    fig3.add_trace(
        go.Scatter(x=Tstr, y=B_clas, name='B classic')
    )
    fig3.add_trace(
        go.Scatter(x=T_state, y=B_virial_state_ref,
                    name='B virial state reference')
        )


    fig1.update_xaxes(range=[0, 10])
    fig1.update_yaxes(range=[-3*De_lim, 3*De_lim])

    fig2.update_xaxes(range=[0, 10])
    fig2.update_yaxes(range=[-3*De_lim, 5*De_lim])



    fig1.update_layout(height=600, width=800,
                    title_text="Graph of Leading Configurations per Distance", xaxis_title=r"$r[A]$",
                    yaxis_title='U[kcal/mol]'
                        )
    st.plotly_chart(fig1, use_container_width=True)

    fig2.update_layout(height=600, width=800,
                    title_text="Graph of Complex Energies per Distance", xaxis_title=r"$r[A]$",
                    yaxis_title='U[kcal/mol]')
    st.plotly_chart(fig2, use_container_width=True)

    fig3.update_layout(height=600, width=800,
                    title_text="Graph of the Second Virial Coefficient B per Temperature T", yaxis_title='B(T)[cm^3/mol]',
                        xaxis_title="T[K]")
    st.plotly_chart(fig3, use_container_width=True)
    if potential == "both":
        fig4 = go.Figure()

        fig4.add_trace(
            go.Scatter(x=r, y=LC_lst_2[0], name='UH')
        )
        fig4.add_trace(
            go.Scatter(x=r, y=LC_lst_2[1], name='UX')
        )
        fig4.add_trace(
            go.Scatter(x=r, y=LC_lst_2[2], name='UZ')
        )
        fig4.add_trace(
            go.Scatter(x=r, y=LC_lst_2[3], name='UTa')
        )
        fig4.add_trace(
            go.Scatter(x=r, y=LC_lst_2[4], name='UTb')
        )
        fig4.add_trace(
            go.Scatter(x=r, y=LC_lst_2[5], name='UL')
        )


        fig5 = go.Figure()
        fig5.add_trace(
            go.Scatter(x=r, y=harm_esf_lst_2[0], name='UM_000')
        )
        fig5.add_trace(
            go.Scatter(x=r, y=harm_esf_lst_2[1], name='UM_202')
        )
        fig5.add_trace(
            go.Scatter(x=r, y=harm_esf_lst_2[2], name='UM_220')
        )
        fig5.add_trace(
            go.Scatter(x=r, y=harm_esf_lst_2[3], name='UM_022')
        )
        fig5.add_trace(
            go.Scatter(x=r, y=harm_esf_lst_2[4], name='UM_222')
        )
        fig5.add_trace(
            go.Scatter(x=r, y=harm_esf_lst_2[5], name='UM_224')
        )



        fig6 = go.Figure()

        fig6.add_trace(
            go.Scatter(x=Tstr_2, y=B_clas_2, name='B classic rydberg')
        )
        fig6.add_trace(
            go.Scatter(x=Tstr, y=B_clas, name='B classic improved lennard-jonnes')
        )
        fig6.add_trace(
        go.Scatter(x=T_state, y=B_virial_state_ref,
                    name='B virial state reference')
        )


        fig4.update_xaxes(range=[0, 10])
        fig4.update_yaxes(range=[-3*De_lim, 3*De_lim])

        fig5.update_xaxes(range=[0, 10])
        fig5.update_yaxes(range=[-3*De_lim, 5*De_lim])



        fig4.update_layout(height=600, width=800,
                        title_text="Graph of Leading Configurations per Distance", xaxis_title=r"$r[A]$",
                        yaxis_title='U[kcal/mol]'
                            )
        st.plotly_chart(fig4, use_container_width=True)

        fig5.update_layout(height=600, width=800,
                        title_text="Graph of Complex Energies per Distance", xaxis_title=r"$r[A]$",
                        yaxis_title='U[kcal/mol]')
        st.plotly_chart(fig5, use_container_width=True)
        fig6.update_layout(height=600, width=800,
                    title_text="Graph of the Second Virial Coefficient B per Temperature T", yaxis_title='B(T)[cm^3/mol]',
                        xaxis_title="T[K]")
        
        st.plotly_chart(fig6, use_container_width=True)


    st.write('evaluation time {}'.format(time.strftime(
        "%H:%M:%S", time.gmtime(time.time()-start))))
    st.write('number of temperatures:', len(Tstr))

    st.success('Calculations finished!', icon="✅")



if st.button('Calculate'):
    start = time.time()
    if uploaded_file is not None:
        
        if potential == 'Improved Leonard-Jonnes Potential':
            result_ILJ = ILJ()
            result_ref = ref_making(result_ILJ.Reqs)

            graph_gen(result_ILJ.LC_lst, result_ILJ.harm_esf_lst, result_ref[0], result_ref[1],
                       result_ILJ.B_clas, result_ILJ.Tstr, result_ILJ.De_lim)
    
            logging.warning('----------------------')
        elif potential == 'Rydberg Potential':
            result_ryd = ryd()
            result_ref = ref_making(result_ryd.Reqs)
            graph_gen(result_ryd.LC_lst, result_ryd.harm_esf_lst, result_ref[0], result_ref[1],
                       result_ryd.B_clas, result_ryd.Tstr, result_ryd.De_lim)
        elif potential == "both":

    
            result_ILJ = ILJ()
            result_ref = ref_making(result_ILJ.Reqs)

            result_ryd = ryd()
            result_ref = ref_making(result_ryd.Reqs)



            result_graph = graph_gen(result_ILJ.LC_lst, result_ILJ.harm_esf_lst, result_ref[0], result_ref[1],
                       result_ILJ.B_clas, result_ILJ.Tstr, result_ILJ.De_lim, result_ryd.LC_lst, result_ryd.harm_esf_lst,
                       result_ryd.B_clas, result_ryd.Tstr, result_ryd.De_lim)
            


       
    else:
        st.info('☝️ Upload a .dat file')
