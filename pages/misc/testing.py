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

logging.basicConfig(filename='log_information.log', encoding='utf-8', level=logging.DEBUG)

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


class ILJ:
    ...

class Data_manage:
    def __init__(self, uploaded_file, st):
        self.uploaded_file = uploaded_file
        self.st = st  # You'll need to provide the streamlit object reference here
        self.T = None
        self.data = None
    def process_data(self):

        data = pd.read_csv(self.uploaded_file, sep="\s+", header=None)
        data.columns = ["alpha", "beta", "mp", "De", "Req"]
        self.st.subheader('DataFrame')
        self.st.write(data)
        self.data = data

        self.h_alpha = self.data.loc[0, 'alpha']
        self.h_beta = self.data.loc[0, 'beta']
        self.h_mp = self.data.loc[0, 'mp']
        self.h_De = self.data.loc[0, 'De'] 
        self.h_Req = self.data.loc[0, 'Req']
            

        self.x_alpha = self.data.loc[1, 'alpha']
        self.x_beta = self.data.loc[1, 'beta']
        self.x_mp = self.data.loc[1, 'mp']
        self.x_De = self.data.loc[1, 'De'] 
        self.x_Req = self.data.loc[1, 'Req']
            

        self.z_alpha = self.data.loc[2, 'alpha']
        self.z_beta = self.data.loc[2, 'beta']
        self.z_mp = self.data.loc[2, 'mp']
        self.z_De = self.data.loc[2, 'De'] 
        self.z_Req = self.data.loc[2, 'Req']
            

        self.ta_alpha = self.data.loc[3, 'alpha']
        self.ta_beta = self.data.loc[3, 'beta']
        self.ta_mp = self.data.loc[3, 'mp']
        self.ta_De = self.data.loc[3, 'De'] 
        self.ta_Req = self.data.loc[3, 'Req']
        

        self.tb_alpha = self.data.loc[4, 'alpha']
        self.tb_beta = self.data.loc[4, 'beta']
        self.tb_mp = self.data.loc[4, 'mp']
        self.tb_De = self.data.loc[4, 'De'] 
        self.tb_Req = self.data.loc[4, 'Req']
            

        self.l_alpha = self.data.loc[5, 'alpha']
        self.l_beta = self.data.loc[5, 'beta']
        self.l_mp = self.data.loc[5, 'mp']
        self.l_De = self.data.loc[5, 'De'] 
        self.l_Req = self.data.loc[5, 'Req']
    

        # Extracting data


class LC:
    def __init__(self, imported):
        self.imported = imported
        
        
    def UH(self, r):
        n = self.imported.h_beta + self.imported.h_alpha * (r / self.imported.h_Req) ** 2
        return self.imported.h_De * ((self.imported.h_mp / (n - self.imported.h_mp)) * (self.imported.h_Req / r) ** n - (n / (n - self.imported.h_mp)) * (self.imported.h_Req / r) ** self.imported.h_mp)

    def UX(self, r):
        n = self.imported.x_beta + self.imported.x_alpha * (r / self.imported.x_Req) ** 2
        return self.imported.x_De * ((self.imported.x_mp / (n - self.imported.x_mp)) * (self.imported.x_Req / r) ** n - (n / (n - self.imported.x_mp)) * (self.imported.x_Req / r) ** self.imported.x_mp)

    def UZ(self, r):
        n = self.imported.z_beta + self.imported.z_alpha * (r / self.imported.z_Req) ** 2
        return self.imported.z_De * ((self.imported.z_mp / (n - self.imported.z_mp)) * (self.imported.z_Req / r) ** n - (n / (n - self.imported.z_mp)) * (self.imported.z_Req / r) ** self.imported.z_mp)

    def UTa(self, r):
        n = self.imported.ta_beta + self.imported.ta_alpha * (r / self.imported.ta_Req) ** 2
        return self.imported.ta_De * ((self.imported.ta_mp / (n - self.imported.ta_mp)) * (self.imported.ta_Req / r) ** n - (n / (n - self.imported.ta_mp)) * (self.imported.ta_Req / r) ** self.imported.ta_mp)

    def UTb(self, r):
        n = self.imported.tb_beta + self.imported.tb_alpha * (r / self.imported.tb_Req) ** 2
        return self.imported.tb_De * ((self.imported.tb_mp / (n - self.imported.tb_mp)) * (self.imported.tb_Req / r) ** n - (n / (n - self.imported.tb_mp)) * (self.imported.tb_Req / r) ** self.imported.tb_mp)

    def UL(self, r):
        n = self.imported.l_beta + self.imported.l_alpha * (r / self.imported.l_Req) ** 2
        return self.imported.l_De * ((self.imported.l_mp / (n - self.imported.l_mp)) * (self.imported.l_Req / r) ** n - (n / (n - self.imported.l_mp)) * (self.imported.l_Req / r) ** self.imported.l_mp)

class Momentum:
    def __init__(self, lc_instance):
        self.lc_instance = lc_instance
        st.write('momentum criada')


        
    def define(self, r):
        self.UH = self.lc_instance.UH(r)
        self.UL = self.lc_instance.UL(r)
        self.UX = self.lc_instance.UX(r)
        self.UZ = self.lc_instance.UZ(r)
        self.UTa = self.lc_instance.UTa(r)
        self.UTb = self.lc_instance.UTb(r)
    
    def UM_000(self, r):
        #st.write('um000 momentum good')
        return (2 * self.UH + self.UL + 2 *
                (self.UTa + self.UTb + self.UX)) / 9

    def UM_202(self, r):
        return 2 * (self.UH - self.UL + self.UTa -
                    2 * self.UTb + self.UX) / (9 * (5**(1 / 2)))

    def UM_022(self, r):
        return 2 * (self.UH - self.UL - 2 * self.UTa +
                    self.UTb + self.UX) / (9 * (5**(1 / 2)))

    def UM_220(self, r):
        return 2 * (4 * self.UH - self.UL - 5 * (self.UTa +
                    self.UTb + self.UX) + 12 * self.UZ) / (45 * (5**(1 / 2)))

    def UM_222(self, r):
        return ((2 / 7)**(1 / 2)) * (13 * self.UH - self.UL + 7 *
                                    (self.UTa + self.UTb - 2 * self.UX) - 12 * self.UZ) / 45

    def UM_224(self, r):
        return ((2 / 35)**(1 / 2) * 8 * self.UH + self.UL + 2 * self.UZ) / 15
    
    '''
    def call_all(self, r):
        UM_000 = self.momentum_instance.UM_000(r)
        UM_202 = self.momentum_instance.UM_202(r)
        UM_022 = self.momentum_instance.UM_022(r)
        UM_220 = self.momentum_instance.UM_220(r)
        UM_222 = self.momentum_instance.UM_222(r)
        UM_224 = self.momentum_instance.UM_224(r)
        return UM_000, UM_202, UM_022
            UM_000, UM_202 = momentum.momentos()

    '''    

class UM_FINAL:
    def __init__(self, momentum_instance):
        self.momentum_instance = momentum_instance
        st.write('UM_FINAL criada')

    def calculate(self, r, th_a, th_b, phi):
        #st.write('calculate um final good')
        self.momentum_instance.define(r)
        UM_000 = self.momentum_instance.UM_000(r)
        UM_202 = self.momentum_instance.UM_202(r)
        UM_022 = self.momentum_instance.UM_022(r)
        UM_220 = self.momentum_instance.UM_220(r)
        UM_222 = self.momentum_instance.UM_222(r)
        UM_224 = self.momentum_instance.UM_224(r)

        UM_FINAL = (
            UM_000 +
            (5 ** (1 / 2)) / 4 * (3 * (np.cos(2 * th_a)) + 1) * UM_202 +
            (5 ** (1 / 2)) / 4 * (3 * (np.cos(2 * th_b)) + 1) * UM_022 +
            (5 ** (1 / 2)) / 16 * UM_220 * (
                (3 * (np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) +
                12 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) +
                3 * (1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi)
            ) -
            (14 ** (1 / 2)) * 5 / 112 * UM_222 * (
                (3 * np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) +
            6 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) -
            3 * (1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi) +
            ((3 * (70) ** (1 / 2)) / 112) * UM_224 * (
                (3 * (np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) -
                8 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) +
                ((1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi)) / 2
            )
        )
        return UM_FINAL

class TemperatureSetter:
    def set_temperature(self, T):
        st.write('t good')
        self.T = T

    

class Integrand:
    def __init__(self, um_final_instance, temperature_instance):
        self.um_final_instance = um_final_instance
        self.temperature_instance = temperature_instance
        st.write('integrand criada')

    def integrand_vegas(self, x):
        self.T = self.temperature_instance.T
        #st.write('integrand good')
        r, th_a, th_b, phi = x
        # Pass LC as an argument to UM_FINAL
        F = self.um_final_instance.calculate(r, th_a, th_b, phi) * constants.rk / self.T
        if F < -1:
            F = -1

        ff = constants.N_A / 4 * np.sin(th_a) * np.sin(th_b) * (r ** 2) * (1 - np.exp(-F))

        return ff
class Mean_req:
    def __init__(self, data):
        self.data = data
        st.write('mean_req criada')

    def calculate_mean(self):
        Reqs = [self.data.loc[0, 'Req'], self.data.loc[1, 'Req'], self.data.loc[2, 'Req'], 
                self.data.loc[3, 'Req'], self.data.loc[0, 'Req'], self.data.loc[0, 'Req']]

        return mean(Reqs)

class Calculate_virial:
    def __init__(self, integrand_instance, mean_req_instance, temperature_instance):
        self.integrand_instance = integrand_instance
        self.mean_req_instance = mean_req_instance
        self.temperature_instance = temperature_instance
        st.write('calculate_virial criada')
    
    def compute_results(self, step):

        #st.write('outisde look good')
        # Integrator and definition of integration limits
        integ = vegas.Integrator(
        [[0.2 * self.mean_req_instance.calculate_mean(), 2 * self.mean_req_instance.calculate_mean()],
        [0, np.pi], [0, np.pi], [0, 2 * np.pi]])
        B_clas = []
        Tstr = []
        B_main = []

        for temp in range(50, 1000, step):
            logging.info(f"Starting integration for temperature {temp}")
            self.temperature_instance.set_temperature(temp)
            result = integ(self.integrand_instance.integrand_vegas, nitn=10, neval=10000)
            B_clas.append(result.mean)
            Tstr.append(temp)
            B_main.append(result.mean)
            logging.info(f"Integration for temperature {temp} completed")

    
        return B_clas, Tstr, B_main
            

if st.button('Calculate'):
    if uploaded_file is not None:
        
        if potential == 'Improved Leonard-Jonnes Potential':
            imported = Data_manage(uploaded_file, st)
            imported.process_data()
            if imported.data is not None:  # Check if data has been processed
                start = time.time()
                lc_instance = LC(imported)
                stop_LC = time.time()
                st.write('LC', stop_LC - start)


                momentum_instance = Momentum(lc_instance)
                stop_momentum = time.time()
                st.write('momentum', stop_momentum - start)

                um_final_instance = UM_FINAL(momentum_instance)
                stop_umfinal = time.time()
                st.write('um_final', stop_umfinal - start)

                temperature_instance = TemperatureSetter()

                integrand_instance = Integrand(um_final_instance, temperature_instance)
                stop_integrand = time.time()
                st.write('integrand', stop_integrand - start)


                mean_req_instance = Mean_req(imported.data)
                st.write(mean_req_instance.calculate_mean())

                virial_calculation_instance = Calculate_virial(integrand_instance, mean_req_instance,temperature_instance)
                stop_virial = time.time()
                st.write('virial_calc', stop_virial - start)

                B_values, T_values, B_main_values = virial_calculation_instance.compute_results(step)
                stop_calc_instance = time.time()
                st.write('LC', stop_calc_instance - start)


                st.write("B values:", B_values)
                st.write("T values:", T_values)
                st.write("Main B values:", B_main_values)
        elif potential == 'Rydberg Potential':
            pass
        elif potential == "both":

            pass


       
    else:
        st.info('☝️ Upload a .dat file')
