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
from constant import constants_fixed, constants_cm, constants_ev, constants_mev, constants_kcalmol

###############################

logging.basicConfig(filename='log_information.log', encoding='utf-8', level=logging.DEBUG)

st.title('Second Virial Coefficient Calculator for A2B2 Molecule')

uploaded_file = st.file_uploader("upload a file with mev as the unit of energy")



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
st.write('note that if you have used the regression page, your data has been converted to meV already!')



rk = constants_mev.rk


class ILJ:
    ...

class Data_manage:
    def __init__(self, uploaded_file, st):
        self.uploaded_file = uploaded_file
        self.st = st  # You'll need to provide the streamlit object reference here
        self.T = None
        self.data_ryd = None
        self.data_ilj = None
    def process_data_ilj(self):

        data = pd.read_csv(self.uploaded_file, sep="\s+", header=None)
        data.columns = ["alpha", "beta", "mp", "De", "Req"]
        self.st.subheader('DataFrame')
        self.st.write(data)
        self.data_ilj = data

        self.h_alpha = self.data_ilj.loc[0, 'alpha']
        self.h_beta = self.data_ilj.loc[0, 'beta']
        self.h_mp = self.data_ilj.loc[0, 'mp']
        self.h_De = self.data_ilj.loc[0, 'De'] 
        self.h_Req = self.data_ilj.loc[0, 'Req']
            
        self.x_alpha = self.data_ilj.loc[1, 'alpha']
        self.x_beta = self.data_ilj.loc[1, 'beta']
        self.x_mp = self.data_ilj.loc[1, 'mp']
        self.x_De = self.data_ilj.loc[1, 'De'] 
        self.x_Req = self.data_ilj.loc[1, 'Req']
            
        self.z_alpha = self.data_ilj.loc[2, 'alpha']
        self.z_beta = self.data_ilj.loc[2, 'beta']
        self.z_mp = self.data_ilj.loc[2, 'mp']
        self.z_De = self.data_ilj.loc[2, 'De'] 
        self.z_Req = self.data_ilj.loc[2, 'Req']
            
        self.ta_alpha = self.data_ilj.loc[3, 'alpha']
        self.ta_beta = self.data_ilj.loc[3, 'beta']
        self.ta_mp = self.data_ilj.loc[3, 'mp']
        self.ta_De = self.data_ilj.loc[3, 'De'] 
        self.ta_Req = self.data_ilj.loc[3, 'Req']
        
        self.tb_alpha = self.data_ilj.loc[4, 'alpha']
        self.tb_beta = self.data_ilj.loc[4, 'beta']
        self.tb_mp = self.data_ilj.loc[4, 'mp']
        self.tb_De = self.data_ilj.loc[4, 'De'] 
        self.tb_Req = self.data_ilj.loc[4, 'Req']
            
        self.l_alpha = self.data_ilj.loc[5, 'alpha']
        self.l_beta = self.data_ilj.loc[5, 'beta']
        self.l_mp = self.data_ilj.loc[5, 'mp']
        self.l_De = self.data_ilj.loc[5, 'De'] 
        self.l_Req = self.data_ilj.loc[5, 'Req']
    

        # Extracting data
    def process_data_ryd(self):
        if potential == "both":
        
            data = pd.read_csv(uploaded_file_2, sep="\s+", header=None)
        else:
            data = pd.read_csv(uploaded_file, sep="\s+", header=None)
        
        data.columns = ["a1", "a2", "a3", "a4", 'a5', 'De', 'Req', 'Eref']
        st.subheader('DataFrame')
        st.write(data)

        self.data_ryd = data
        self.h_a1 = self.data_ryd.loc[0, 'a1']
        self.h_a2 = self.data_ryd.loc[0, 'a2']
        self.h_a3 = self.data_ryd.loc[0, 'a3']
        self.h_a4 = self.data_ryd.loc[0, 'a4']
        self.h_a5 = self.data_ryd.loc[0, 'a5'] 
        self.h_De = self.data_ryd.loc[0, 'De'] 
        self.h_Req = self.data_ryd.loc[0, 'Req']
        self.h_Eref = self.data_ryd.loc[0, 'Eref']

        self.x_a1 = self.data_ryd.loc[1, 'a1']
        self.x_a2 = self.data_ryd.loc[1, 'a2']
        self.x_a3 = self.data_ryd.loc[1, 'a3']
        self.x_a4 = self.data_ryd.loc[1, 'a4']
        self.x_a5 = self.data_ryd.loc[1, 'a5']
        self.x_De = self.data_ryd.loc[1, 'De']  # eV
        self.x_Req = self.data_ryd.loc[1, 'Req']
        self.x_Eref = self.data_ryd.loc[1, 'Eref']

        self.z_a1 = self.data_ryd.loc[2, 'a1']
        self.z_a2 = self.data_ryd.loc[2, 'a2']
        self.z_a3 = self.data_ryd.loc[2, 'a3']
        self.z_a4 = self.data_ryd.loc[2, 'a4']
        self.z_a5 = self.data_ryd.loc[2, 'a5']
        self.z_De = self.data_ryd.loc[2, 'De']
        self.z_Req = self.data_ryd.loc[2, 'Req']
        self.z_Eref =self.data_ryd.loc[2, 'Eref']

        self.ta_a1 = self.data_ryd.loc[3, 'a1']
        self.ta_a2 = self.data_ryd.loc[3, 'a2']
        self.ta_a3 = self.data_ryd.loc[3, 'a3']
        self.ta_a4 = self.data_ryd.loc[3, 'a4']
        self.ta_a5 = self.data_ryd.loc[3, 'a5']
        self.ta_De = self.data_ryd.loc[3, 'De']
        self.ta_Req = self.data_ryd.loc[3, 'Req']
        self.ta_Eref = self.data_ryd.loc[3, 'Eref']

        self.tb_a1 = self.data_ryd.loc[4, 'a1']
        self.tb_a2 = self.data_ryd.loc[4, 'a2']
        self.tb_a3 = self.data_ryd.loc[4, 'a3']
        self.tb_a4 = self.data_ryd.loc[4, 'a4']
        self.tb_a5 = self.data_ryd.loc[4, 'a5']
        self.tb_De = self.data_ryd.loc[4, 'De']
        self.tb_Req = self.data_ryd.loc[4, 'Req']
        self.tb_Eref = self.data_ryd.loc[4, 'Eref']

        self.l_a1 = self.data_ryd.loc[5, 'a1']
        self.l_a2 = self.data_ryd.loc[5, 'a2']
        self.l_a3 = self.data_ryd.loc[5, 'a3']
        self.l_a4 = self.data_ryd.loc[5, 'a4']
        self.l_a5 = self.data_ryd.loc[5, 'a5']
        self.l_De = self.data_ryd.loc[5, 'De']
        self.l_Req = self.data_ryd.loc[5, 'Req']
        self.l_Eref = self.data_ryd.loc[5, 'Eref']


class LCs_ilj:
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



class LCs_ryd:
    def __init__(self, imported):
        self.imported = imported

    def UH(self, r):
        y = (r - self.imported.h_Req) / self.imported.h_Req
        return (-self.imported.h_De * (1 + self.imported.h_a1 * y + self.imported.h_a2 * y**2 + self.imported.h_a3 * y**3 + self.imported.h_a4 * y**4 + self.imported.h_a5 * y**5) * np.exp(-self.imported.h_a1 * y) + self.imported.h_Eref)

    def UX(self, r):
        y = (r - self.imported.x_Req) / self.imported.x_Req
        return (-self.imported.x_De * (1 + self.imported.x_a1 * y + self.imported.x_a2 * y**2 + self.imported.x_a3 * y**3 + self.imported.x_a4 * y**4 + self.imported.x_a5 * y**5) * np.exp(-self.imported.x_a1 * y) + self.imported.x_Eref)
        

    def UZ(self, r):
        y = (r - self.imported.z_Req) / self.imported.z_Req
        return (-self.imported.z_De * (1 + self.imported.z_a1 * y + self.imported.z_a2 * y**2 + self.imported.z_a3 * y**3 + self.imported.z_a4 * y**4 + self.imported.z_a5 * y**5) * np.exp(-self.imported.z_a1 * y) + self.imported.z_Eref)

    def UTa(self, r):
        y = (r - self.imported.ta_Req) / self.imported.ta_Req
        return (-self.imported.ta_De * (1 + self.imported.ta_a1 * y + self.imported.ta_a2 * y**2 + self.imported.ta_a3 * y**3 + self.imported.ta_a4 * y**4 + self.imported.ta_a5 * y**5) * np.exp(-self.imported.ta_a1 * y) + self.imported.ta_Eref)
        

    def UTb(self, r):
        y = (r - self.imported.tb_Req) / self.imported.tb_Req
        return (-self.imported.tb_De * (1 + self.imported.tb_a1 * y + self.imported.tb_a2 * y**2 + self.imported.tb_a3 * y**3 + self.imported.tb_a4 * y**4 + self.imported.tb_a5 * y**5) * np.exp(-self.imported.tb_a1 * y) + self.imported.tb_Eref)

    def UL(self, r):
        y = (r - self.imported.l_Req) / self.imported.l_Req
        return (-self.imported.l_De * (1 + self.imported.l_a1 * y + self.imported.l_a2 * y**2 + self.imported.l_a3 * y**3 + self.imported.l_a4 * y**4 + self.imported.l_a5 * y**5) * np.exp(-self.imported.l_a1 * y) + self.imported.l_Eref)
        

    

class momentum_ilj(LCs_ilj):
    def __init__(self):
        super().__init__(imported)

    def UM_000(self, r):
        return (2 * self.UH(r) + self.UL(r) + 2 * (self.UTa(r) + self.UTb(r) + self.UX(r))) / 9

    def UM_202(self, r):
        return 2 * (self.UH(r) - self.UL(r) + self.UTa(r) -
                    2 * self.UTb(r) + self.UX(r)) / (9 * (5**(1 / 2)))

    def UM_022(self, r):
        return 2 * (self.UH(r) - self.UL(r) - 2 * self.UTa(r) +
                    self.UTb(r) + self.UX(r)) / (9 * (5**(1 / 2)))

    def UM_220(self, r):
        return 2 * (4 * self.UH(r) - self.UL(r) - 5 * (self.UTa(r) +
                    self.UTb(r) + self.UX(r)) + 12 * self.UZ(r)) / (45 * (5**(1 / 2)))

    def UM_222(self, r):
        return ((2 / 7)**(1 / 2)) * (13 * self.UH(r) - self.UL(r) + 7 *
                                    (self.UTa(r) + self.UTb(r) - 2 * self.UX(r)) - 12 * self.UZ(r)) / 45

    def UM_224(self, r):
        return ((2 / 35)**(1 / 2) * 8 * self.UH(r) + self.UL(r) + 2 * self.UZ(r)) / 15
    
class momentum_ryd(LCs_ryd):
    def __init__(self):
        super().__init__(imported)

    def UM_000(self, r):
        return (2 * self.UH(r) + self.UL(r) + 2 * (self.UTa(r) + self.UTb(r) + self.UX(r))) / 9

    def UM_202(self, r):
        return 2 * (self.UH(r) - self.UL(r) + self.UTa(r) -
                    2 * self.UTb(r) + self.UX(r)) / (9 * (5**(1 / 2)))

    def UM_022(self, r):
        return 2 * (self.UH(r) - self.UL(r) - 2 * self.UTa(r) +
                    self.UTb(r) + self.UX(r)) / (9 * (5**(1 / 2)))

    def UM_220(self, r):
        return 2 * (4 * self.UH(r) - self.UL(r) - 5 * (self.UTa(r) +
                    self.UTb(r) + self.UX(r)) + 12 * self.UZ(r)) / (45 * (5**(1 / 2)))

    def UM_222(self, r):
        return ((2 / 7)**(1 / 2)) * (13 * self.UH(r) - self.UL(r) + 7 *
                                    (self.UTa(r) + self.UTb(r) - 2 * self.UX(r)) - 12 * self.UZ(r)) / 45

    def UM_224(self, r):
        return ((2 / 35)**(1 / 2) * 8 * self.UH(r) + self.UL(r) + 2 * self.UZ(r)) / 15


class derivatives_ilj(LCs_ilj):
    def __init__(self):
        super().__init__(imported)


    def d1_H(self, r):
        d1_H = self.derivative(
            self.UH, r, method='central', p=1e-8)
        return d1_H

    def d1_L(self, r):
        d1_L = self.derivative(
            self.UL, r, method='central', p=1e-8)
        return d1_L

    def d1_Ta(self, r):
        d1_Ta = self.derivative(
            self.UTa, r, method='central', p=1e-8)
        return d1_Ta

    def d1_Tb(self, r):
        d1_Tb = self.derivative(
            self.UTb, r, method='central', p=1e-8)
        return d1_Tb

    def d1_X(self, r):
        d1_X = self.derivative(
            self.UX, r, method='central', p=1e-8)
        return d1_X

    def d1_Z(self, r):
        d1_Z = self.derivative(
            self.UZ, r, method='central', p=1e-8)
        return d1_Z

    def d2_H(self, r):
        d2_H = self.derivative2(
            self.UH, r, method='central', p=1e-8)
        return d2_H

    def d2_L(self, r):
        d2_L = self.derivative2(
            self.UL, r, method='central', p=1e-8)
        return d2_L

    def d2_Ta(self, r):
        d2_Ta = self.derivative2(
            self.UTa, r, method='central', p=1e-8)
        return d2_Ta

    def d2_Tb(self, r):
        d2_Tb = self.derivative2(
            self.UTb, r, method='central', p=1e-8)
        return d2_Tb

    def d2_X(self, r):
        d2_X = self.derivative2(
            self.UX, r, method='central', p=1e-8)
        return d2_X

    def d2_Z(self, r):
        d2_Z = self.derivative2(
            self.UZ, r, method='central', p=1e-8)
        return d2_Z
    
    def d1r(self, r, th_a, th_b, phi):
        d1r_1 = ((-1/18)*(1+3*np.cos(2*th_a))*(self.d1_H(r) - self.d1_L(r)+self.d1_Ta(r)-2*self.d1_Tb(r)+self.d1_X(r)))
        d1r_2 = (-(1/18)*(1+3*np.cos(2*th_b))*(self.d1_H(r) -self.d1_L(r)-2*self.d1_Ta(r)+self.d1_Tb(r)+self.d1_X(r)))
        d1r_3 = (+(1/9)*(2*self.d1_H(r)+self.d1_L(r)+2 *(self.d1_Ta(r)+self.d1_Tb(r)+self.d1_X(r))))
        d1r_4 = (-(1/504)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))- 3*(1-np.cos(2*th_a)) *(1-np.cos(2*th_b))*np.cos(2*phi)+ 6*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(13*self.d1_H(r)- self.d1_L(r)+7*self.d1_Ta(r)+7*self.d1_Tb(r)-14*self.d1_X(r)-12*self.d1_Z(r)))
        d1r_5 = (+(1/35)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+1/2*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)- 8*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(self.d1_H(r)+self.d1_L(r)-2*self.d1_Z(r)))
        d1r_6 = (+(1/360)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)+ 12*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(4*self.d1_H(r)-self.d1_L(r)-5*self.d1_Ta(r)-5*self.d1_Tb(r)-5*self.d1_X(r)+12*self.d1_Z(r)))
        return d1r_1 + d1r_2 + d1r_3 + d1r_4 + d1r_5 + d1r_6

    def d1th_a(self, r, th_a, th_b, phi):
        d1th_a_1 = (1/3*np.sin(2*th_a) * (self.UH(r) - self.UL(r)+self.UTa(r)-2*self.UTb(r)+self.UX(r)))
        d1th_a_2 = (-(1/504)*(-6*(1+3*np.cos(2*th_b)) * np.sin(2*th_a) - 6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)+ 12*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b))* (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r)))
        d1th_a_3 = (+(1/35)*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + (1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a) - 16*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (self.UH(r)+self.UL(r)-2*self.UZ(r)))
        d1th_a_4 = (+(1/360)*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + 6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a) + 24*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r)))
        return d1th_a_1 + d1th_a_2 + d1th_a_3 + d1th_a_4

    def d1th_b(self, r, th_a, th_b, phi):
        d1th_b_1 = 1/3*np.sin(2*th_b) * (self.UH(r) -self.UL(r)-2*self.UTa(r)+self.UTb(r)+self.UX(r))
        d1th_b_2 = -1/504*(12*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a) - 6*(1+3*np.cos(2*th_a)) * np.sin(2*th_b) - 6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d1th_b_3 = +1/35*(-16*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a)- 6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b) + (1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(self.UH(r)+self.UL(r)-2*self.UZ(r))
        d1th_b_4 = +1/360*(24*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a) - 6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b)+ 6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d1th_b_1 + d1th_b_2 + d1th_b_3 + d1th_b_4

    def d1phi(self,r, th_a, th_b, phi):
        d1phi_1 = -1/504*(-6*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) + 6*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi)) * (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d1phi_2 = +1/35*(8*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - (1-np.cos(2*th_a))* (1-np.cos(2*th_b)) * np.sin(2*phi)) * (self.UH(r)+self.UL(r)-2*self.UZ(r))
        d1phi_3 = + 1/360*(-12*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - 6*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi))*(4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d1phi_1 + d1phi_2 + d1phi_3

    def d2r(self, r, th_a, th_b, phi):
        d2r_1 = -1/18*(1+3*np.cos(2*th_a)) * (self.d2_H(r) -self.d2_L(r)+self.d2_Ta(r)-2*self.d2_Tb(r)+self.d2_X(r))
        d2r_2 = -1/18*(1+3*np.cos(2*th_b)) * (self.d2_H(r) -self.d2_L(r)-2*self.d2_Ta(r)+self.d2_Tb(r)+self.d2_X(r))
        d2r_3 = +1/9 * (2*self.d2_H(r)+self.d2_L(r)+2 *(self.d2_Ta(r)+self.d2_Tb(r)+self.d2_X(r)))
        d2r_4 = -1/504*(1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b))
        d2r_5 = -3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi) + 6*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b) * (13*self.d2_H(r)-self.d2_L(r)+7*self.d2_Ta(r)+7*self.d2_Tb(r)-14*self.d2_X(r)-12*self.d2_Z(r))
        d2r_6 = +1/35*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 1/2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)- 8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (self.d2_H(r)+self.d2_L(r)-2*self.d2_Z(r))
        d2r_7 = +1/360*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*self.d2_H(r)-self.d2_L(r)-5*self.d2_Ta(r)-5*self.d2_Tb(r)-5*self.d2_X(r)+12*self.d2_Z(r))
        return d2r_1 + d2r_2 + d2r_3 + d2r_4 + d2r_5 + d2r_6 + d2r_7

    def d2th_a(self, r, th_a, th_b, phi):
        d2th_a_1 = 2/3*np.cos(2*th_a) * (self.UH(r) -self.UL(r)+self.UTa(r)-2*self.UTb(r)+self.UX(r))
        d2th_a_2 = + 1/504*(12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d2th_a_3 = + 1/35*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 2*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (self.UH(r)+self.UL(r)-2*self.UZ(r))
        d2th_a_4 = + 1/360*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)- 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d2th_a_1 + d2th_a_2 + d2th_a_3 + d2th_a_4

    def d2th_b(self, r, th_a, th_b, phi):
        d2th_b_1 = 2/3*np.cos(2*th_b) * (self.UH(r) - self.UL(r)-2*self.UTa(r)+self.UTb(r)+self.UX(r))
        d2th_b_2 = + 1/504*(12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 12*(1-np.cos(2*th_a)) * np.cos(2*th_b) * np.cos(2*phi)+ 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d2th_b_3 = + 1/35*(-12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 2*(1-np.cos(2*th_a) * np.cos(2*th_b)) * np.cos(2*phi)+ 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (self.UH(r)+self.UL(r)-2*self.UZ(r))
        d2th_b_4 = + 1/360*(-12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 12*(1-np.cos(2*th_a)) * np.cos(2*th_b) * np.cos(2*phi)- 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d2th_b_1 + d2th_b_2 + d2th_b_3 + d2th_b_4

    def d2phi(self, r, th_a, th_b, phi):
        d2phi_1 = 1/504*(-12*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 6*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d2phi_2 = + 1/35*(-2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (self.UH(r)+self.UL(r)-2*self.UZ(r))
        d2phi_3 = + 1/360*(-12*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)- 12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d2phi_1 + d2phi_2 + d2phi_3
    def derivative(self, g, a, method='central', p=0.01):
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

    def derivative2(self, g, a, method='central', p=0.01):  # derivada de segunda ordem
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
    
    
class derivatives_ryd(LCs_ryd):
    def __init__(self):
        super().__init__(imported)


    def d1_H(self, r):
        d1_H = self.derivative(
            self.UH, r, method='central', p=1e-8)
        return d1_H

    def d1_L(self, r):
        d1_L = self.derivative(
            self.UL, r, method='central', p=1e-8)
        return d1_L

    def d1_Ta(self, r):
        d1_Ta = self.derivative(
            self.UTa, r, method='central', p=1e-8)
        return d1_Ta

    def d1_Tb(self, r):
        d1_Tb = self.derivative(
            self.UTb, r, method='central', p=1e-8)
        return d1_Tb

    def d1_X(self, r):
        d1_X = self.derivative(
            self.UX, r, method='central', p=1e-8)
        return d1_X

    def d1_Z(self, r):
        d1_Z = self.derivative(
            self.UZ, r, method='central', p=1e-8)
        return d1_Z

    def d2_H(self, r):
        d2_H = self.derivative2(
            self.UH, r, method='central', p=1e-8)
        return d2_H

    def d2_L(self, r):
        d2_L = self.derivative2(
            self.UL, r, method='central', p=1e-8)
        return d2_L

    def d2_Ta(self, r):
        d2_Ta = self.derivative2(
            self.UTa, r, method='central', p=1e-8)
        return d2_Ta

    def d2_Tb(self, r):
        d2_Tb = self.derivative2(
            self.UTb, r, method='central', p=1e-8)
        return d2_Tb

    def d2_X(self, r):
        d2_X = self.derivative2(
            self.UX, r, method='central', p=1e-8)
        return d2_X

    def d2_Z(self, r):
        d2_Z = self.derivative2(
            self.UZ, r, method='central', p=1e-8)
        return d2_Z

    
    def d1r(self, r, th_a, th_b, phi):
        d1r_1 = ((-1/18)*(1+3*np.cos(2*th_a))*(self.d1_H(r) - self.d1_L(r)+self.d1_Ta(r)-2*self.d1_Tb(r)+self.d1_X(r)))
        d1r_2 = (-(1/18)*(1+3*np.cos(2*th_b))*(self.d1_H(r) -self.d1_L(r)-2*self.d1_Ta(r)+self.d1_Tb(r)+self.d1_X(r)))
        d1r_3 = (+(1/9)*(2*self.d1_H(r)+self.d1_L(r)+2 *(self.d1_Ta(r)+self.d1_Tb(r)+self.d1_X(r))))
        d1r_4 = (-(1/504)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))- 3*(1-np.cos(2*th_a)) *(1-np.cos(2*th_b))*np.cos(2*phi)+ 6*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(13*self.d1_H(r)- self.d1_L(r)+7*self.d1_Ta(r)+7*self.d1_Tb(r)-14*self.d1_X(r)-12*self.d1_Z(r)))
        d1r_5 = (+(1/35)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+1/2*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)- 8*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(self.d1_H(r)+self.d1_L(r)-2*self.d1_Z(r)))
        d1r_6 = (+(1/360)*((1+3*np.cos(2*th_a))*(1+3*np.cos(2*th_b))+3*(1-np.cos(2*th_a))*(1-np.cos(2*th_b))*np.cos(2*phi)+ 12*np.cos(phi)*np.sin(2*th_a)*np.sin(2*th_b))*(4*self.d1_H(r)-self.d1_L(r)-5*self.d1_Ta(r)-5*self.d1_Tb(r)-5*self.d1_X(r)+12*self.d1_Z(r)))
        return d1r_1 + d1r_2 + d1r_3 + d1r_4 + d1r_5 + d1r_6

    def d1th_a(self, r, th_a, th_b, phi):
        d1th_a_1 = (1/3*np.sin(2*th_a) * (self.UH(r) - self.UL(r)+self.UTa(r)-2*self.UTb(r)+self.UX(r)))
        d1th_a_2 = (-(1/504)*(-6*(1+3*np.cos(2*th_b)) * np.sin(2*th_a) - 6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a)+ 12*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b))* (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r)))
        d1th_a_3 = (+(1/35)*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + (1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a) - 16*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (self.UH(r)+self.UL(r)-2*self.UZ(r)))
        d1th_a_4 = (+(1/360)*(-6 * (1+3*np.cos(2*th_b)) * np.sin(2*th_a) + 6*(1-np.cos(2*th_b)) * np.cos(2*phi) * np.sin(2*th_a) + 24*np.cos(2*th_a) * np.cos(phi) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r)))
        return d1th_a_1 + d1th_a_2 + d1th_a_3 + d1th_a_4

    def d1th_b(self, r, th_a, th_b, phi):
        d1th_b_1 = 1/3*np.sin(2*th_b) * (self.UH(r) -self.UL(r)-2*self.UTa(r)+self.UTb(r)+self.UX(r))
        d1th_b_2 = -1/504*(12*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a) - 6*(1+3*np.cos(2*th_a)) * np.sin(2*th_b) - 6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d1th_b_3 = +1/35*(-16*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a)- 6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b) + (1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b))*(self.UH(r)+self.UL(r)-2*self.UZ(r))
        d1th_b_4 = +1/360*(24*np.cos(2*th_b) * np.cos(phi) * np.sin(2*th_a) - 6 * (1+3*np.cos(2*th_a)) * np.sin(2*th_b)+ 6*(1-np.cos(2*th_a)) * np.cos(2*phi) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d1th_b_1 + d1th_b_2 + d1th_b_3 + d1th_b_4

    def d1phi(self,r, th_a, th_b, phi):
        d1phi_1 = -1/504*(-6*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) + 6*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi)) * (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d1phi_2 = +1/35*(8*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - (1-np.cos(2*th_a))* (1-np.cos(2*th_b)) * np.sin(2*phi)) * (self.UH(r)+self.UL(r)-2*self.UZ(r))
        d1phi_3 = + 1/360*(-12*np.sin(2*th_a) * np.sin(2*th_b) * np.sin(phi) - 6*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.sin(2*phi))*(4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d1phi_1 + d1phi_2 + d1phi_3

    def d2r(self, r, th_a, th_b, phi):
        d2r_1 = -1/18*(1+3*np.cos(2*th_a)) * (self.d2_H(r) -self.d2_L(r)+self.d2_Ta(r)-2*self.d2_Tb(r)+self.d2_X(r))
        d2r_2 = -1/18*(1+3*np.cos(2*th_b)) * (self.d2_H(r) -self.d2_L(r)-2*self.d2_Ta(r)+self.d2_Tb(r)+self.d2_X(r))
        d2r_3 = +1/9 * (2*self.d2_H(r)+self.d2_L(r)+2 *(self.d2_Ta(r)+self.d2_Tb(r)+self.d2_X(r)))
        d2r_4 = -1/504*(1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b))
        d2r_5 = -3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi) + 6*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b) * (13*self.d2_H(r)-self.d2_L(r)+7*self.d2_Ta(r)+7*self.d2_Tb(r)-14*self.d2_X(r)-12*self.d2_Z(r))
        d2r_6 = +1/35*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 1/2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)- 8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (self.d2_H(r)+self.d2_L(r)-2*self.d2_Z(r))
        d2r_7 = +1/360*((1+3*np.cos(2*th_a)) * (1+3*np.cos(2*th_b)) + 3*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*self.d2_H(r)-self.d2_L(r)-5*self.d2_Ta(r)-5*self.d2_Tb(r)-5*self.d2_X(r)+12*self.d2_Z(r))
        return d2r_1 + d2r_2 + d2r_3 + d2r_4 + d2r_5 + d2r_6 + d2r_7

    def d2th_a(self, r, th_a, th_b, phi):
        d2th_a_1 = 2/3*np.cos(2*th_a) * (self.UH(r) -self.UL(r)+self.UTa(r)-2*self.UTb(r)+self.UX(r))
        d2th_a_2 = + 1/504*(12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d2th_a_3 = + 1/35*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 2*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (self.UH(r)+self.UL(r)-2*self.UZ(r))
        d2th_a_4 = + 1/360*(-12*np.cos(2*th_a) * (1+3*np.cos(2*th_b)) + 12*np.cos(2*th_a) * (1-np.cos(2*th_b)) * np.cos(2*phi)- 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d2th_a_1 + d2th_a_2 + d2th_a_3 + d2th_a_4

    def d2th_b(self, r, th_a, th_b, phi):
        d2th_b_1 = 2/3*np.cos(2*th_b) * (self.UH(r) - self.UL(r)-2*self.UTa(r)+self.UTb(r)+self.UX(r))
        d2th_b_2 = + 1/504*(12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 12*(1-np.cos(2*th_a)) * np.cos(2*th_b) * np.cos(2*phi)+ 24*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d2th_b_3 = + 1/35*(-12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 2*(1-np.cos(2*th_a) * np.cos(2*th_b)) * np.cos(2*phi)+ 32*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (self.UH(r)+self.UL(r)-2*self.UZ(r))
        d2th_b_4 = + 1/360*(-12*(1+3*np.cos(2*th_a)) * np.cos(2*th_b) + 12*(1-np.cos(2*th_a)) * np.cos(2*th_b) * np.cos(2*phi)- 48*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d2th_b_1 + d2th_b_2 + d2th_b_3 + d2th_b_4

    def d2phi(self, r, th_a, th_b, phi):
        d2phi_1 = 1/504*(-12*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 6*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (13*self.UH(r)-self.UL(r)+7*self.UTa(r)+7*self.UTb(r)-14*self.UX(r)-12*self.UZ(r))
        d2phi_2 = + 1/35*(-2*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)+ 8*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (self.UH(r)+self.UL(r)-2*self.UZ(r))
        d2phi_3 = + 1/360*(-12*(1-np.cos(2*th_a)) * (1-np.cos(2*th_b)) * np.cos(2*phi)- 12*np.cos(phi) * np.sin(2*th_a) * np.sin(2*th_b)) * (4*self.UH(r)-self.UL(r)-5*self.UTa(r)-5*self.UTb(r)-5*self.UX(r)+12*self.UZ(r))
        return d2phi_1 + d2phi_2 + d2phi_3
    def derivative(self, g, a, method='central', p=0.01):
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

    def derivative2(self, g, a, method='central', p=0.01):  # derivada de segunda ordem
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
    

 

class UM_FINAL:
    def __init__(self, momentum_instance):
        self.momentum_instance = momentum_instance
        #st.write('UM_FINAL criada')

    def calculate(self, r, th_a, th_b, phi):
        #st.write('calculate um final good')
        UM_000 = self.momentum_instance.UM_000(r)
        UM_202 = self.momentum_instance.UM_202(r)
        UM_022 = self.momentum_instance.UM_022(r)
        UM_220 = self.momentum_instance.UM_220(r)
        UM_222 = self.momentum_instance.UM_222(r)
        UM_224 = self.momentum_instance.UM_224(r)

        UM_FINAL_1 = UM_000 + (5 ** (1 / 2)) / 4 * (3 * (np.cos(2 * th_a)) + 1) * UM_202 
        UM_FINAL_2 = +(5 ** (1 / 2)) / 4 * (3 * (np.cos(2 * th_b)) + 1) * UM_022 
        UM_FINAL_3 = +(5 ** (1 / 2)) / 16 * UM_220 * ((3 * (np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) +12 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) +3 * (1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi))
        UM_FINAL_4 = -(14 ** (1 / 2)) * 5 / 112 * UM_222 * ((3 * np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) +6 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) -3 * (1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi) 
        UM_FINAL_5 = +((3 * (70) ** (1 / 2)) / 112) * UM_224 * ((3 * (np.cos(2 * th_a)) + 1) * (3 * (np.cos(2 * th_b)) + 1) -8 * np.sin(2 * th_a) * np.sin(2 * th_b) * np.cos(phi) +((1 - np.cos(2 * th_a)) * (1 - np.cos(2 * th_b)) * np.cos(2 * phi)) / 2)
        
        return UM_FINAL_1 + UM_FINAL_2 + UM_FINAL_3 + UM_FINAL_4 + UM_FINAL_5
    
    

    
class RefMaking:
    def __init__(self, gas, media_Reqs, temperature_instance):
        self.gas = gas
        self.media_Reqs =   media_Reqs
        self.temperature_instance = temperature_instance
        self.massaC = 12.0
        self.massaO = 15.99491
        self.massaN = 14.00307
        self.massaH = 1.00783
        self.massaF = 18.99840
        self.massaBr = 79.904
        self.massaCl = 34.96885
        self.mi = None
        self.I1 = None
        self.I2 = None
        self.mi1 = None
        self.mi2 = None
        self.B_virial_state_ref = []
        self.T_state = []



    def BH2_ref(self,T):  # from 60K to 500K
        return 1.7472*10 - (1.2926 * 10**2)/(T) - (2.6988 * 10**5)/(T)**2 + (8.0282 * 10**6)/(T)**3

    def BCl2_ref(self,T):  # from  300K  to 1070?
        return 1.3171*10 + 5.1055*10**3/T - 2.9404*10**7/T**2 - 5.0383*10**8/T**3

    def BF2_ref(self,T):  # from  85K  to 295?
        return 3.3609*10 - 1.0625*10**4/T - 6.0780*10**5/T**2 - 2.2759*10**7/T**3
        

    

    def gen_ref(self):
    
        if gas == 'H2F2':
            #st.write('h2f2')
            m1 = self.massaH
            m2 = self.massaH
            m3 = self.massaF
            m4 = self.massaF

            Req1 = 0.744013
            Req2 = 1.415962  # Cl  2.007880 # F = 1.415962 #usar dados CBS da tabela
            T_ref = [80, 90, 100, 110, 120, 160, 200, 250, 300]
            B_ref = [-50.13, -50.13, -33.91, -27.86,  -
                        22.83,  -9.18,  -1.30,  4.70, 8.49]

            self.mi = (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
            #mi1 = (m1 + m2) / 4 * Req1**2 / 6.02252
            #mi2 = (m3 + m4) / 4 * Req2**2 / 6.02252

            self.r1 = (m2 * Req1) / (m1+m2) #(m2 * self.media_Reqs) / (m1 - m2)
            self.r2 = Req1 -  self.r1
            self.r3 = (m4 * Req2) / (m3+m4) #(m4 * self.media_Reqs) / (m3 - m4)
            self.r4 = Req2 -  self.r3

            self.I1 = m1 *  self.r1**2 + m2 *  self.r2**2
            self.I2 =  m3 *  self.r3**2 + m4 *  self.r4**2


            self.mi1 = (m1 * m2)/(m1 + m2)
            self.mi2 = (m3 * m4)/(m3 + m4)                


            self.Ma = m1 + m2
            self.Mb = m3 + m4
            self.w1 = -0.216
            self.w2 = 0.053
            self.rho_1 = 0.026667
            self.rho_2 = 0.02715
            self.zc1 = 0.303
            self.zc2 = 0.287
            self.x1 = 0.5
            self.x2 = 0.5

            def B_A2_ref(T):
                return self.BF2_ref(T)

            def B_B2_ref(T):
                return self.BH2_ref(T)

            self.B_A2_ref = B_A2_ref
            self.B_B2_ref = B_B2_ref
        
        elif gas == 'F2F2':
            #st.write('f2f2')
            m1 = self.massaF
            m2 = self.massaF
            m3 = self.massaF
            m4 = self.massaF

            Req1 = 1.415962
            Req2 = 1.415962  # Cl  2.007880 # F = 1.415962 #usar dados CBS da tabela

            T_ref = [80, 90, 100, 110, 120, 160, 200, 250, 300]  # errado
            B_ref = [-50.13, -50.13, -33.91, -27.86,  -22.83,  -
                        9.18,  -1.30,  4.70, 8.49]  # errado, sao do h2f2
            
            self.mi =  (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
            #mi1 = (m1 + m2) / 4 * Req1**2 / 6.02252
            #mi2 = (m3 + m4) / 4 * Req2**2 / 6.02252
            
            
            
            self.r1 = (m2 * Req1) / (m1+m2) #(m2 * self.media_Reqs) / (m1 - m2)
            self.r2 = Req1 - self.r1
            self.r3 = (m4 * Req2) / (m3+m4) #(m4 * self.media_Reqs) / (m3 - m4)
            self.r4 = Req2 - self.r3


            self.I1 = m1 * self.r1**2 + m2 *  self.r2**2
            self.I2 =  m3 *  self.r3**2 + m4 *  self.r4**2


            self.mi1 = (m1 * m2)/(m1 + m2)
            self.mi2 = (m3 * m4)/(m3 + m4)


            self.Ma = m1 + m2
            self.Mb = m3 + m4
            self.w1 = 0.053
            self.w2 = 0.053
            self.rho_1 = 0.02715
            self.rho_2 = 0.02715
            self.zc1 = 0.287
            self.zc2 = 0.287
            self.x1 = 0.5
            self.x2 = 0.5

            def B_A2_ref(T):
                return self.BF2_ref(T)

            def B_B2_ref(T):
                return self.BF2_ref(T)
            
            self.B_A2_ref = B_A2_ref
            self.B_B2_ref = B_B2_ref
        elif gas == 'H2H2':
            #st.write('h2h2')

            m1 = self.massaH
            m2 = self.massaH
            m3 = self.massaH
            m4 = self.massaH
            Req1 = 0.744013
            Req2 = 0.744013

            self.mi =  (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
            #mi1 = (m1 + m2) / 4 * Req1**2 / 6.02252
            #mi2 = (m3 + m4) / 4 * Req2**2 / 6.02252

            self.r1 = (m2 * Req1) / (m1+m2) #(m2 * self.media_Reqs) / (m1 - m2)
            self.r2 = Req1 - self.r1
            self.r3 = (m4 * Req2) / (m3+m4) #(m4 * self.media_Reqs) / (m3 - m4)
            self.r4 = Req2 - self.r3



            self.I1 = m1 * self.r1**2 + m2 * self.r2**2
            self.I2 =  m3 * self.r3**2 + m4 * self.r4**2


            self.mi1 = (m1 * m2)/(m1 + m2)
            self.mi2 = (m3 * m4)/(m3 + m4)

            self.Ma = m1 + m2
            self.Mb = m3 + m4
            self.w1 = -0.216
            self.w2 = -0.216
            self.rho_1 = 0.026667
            self.rho_2 = 0.026667
            self.zc1 = 0.303
            self.zc2 = 0.303
            self.x1 = 0.5
            self.x2 = 0.5

            def B_A2_ref(T):
                return self.BH2_ref(T)

            def B_B2_ref(T):
                return self.BH2_ref(T)
            self.B_A2_ref = B_A2_ref
            self.B_B2_ref = B_B2_ref

        elif gas == 'H2Br2':
            pass
            '''
            #st.write('h2Br2')
            m1 = self.massaH
            m2 = self.massaH
            m3 = self.massaBr
            m4 = self.massaBr
            B_ref = []
            T_ref = []
            Req1 = 0
            Req2 = 0

            mi = (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
            #mi1 = (m1 + m2) / 4 * Req1**2 / 6.02252
            #mi2 = (m3 + m4) / 4 * Req2**2 / 6.02252

            r1 = (m2 * self.media_Reqs) / (m1 - m2)
            r2 = self.media_Reqs - r1
            r3 = (m4 * self.media_Reqs) / (m3 - m4)
            r4 = self.media_Reqs - r3



            self.I1 = m1 * r1**2 + m2 * r2**2
            self.I2 =  m3 * r3**2 + m4 * r4**2


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
            self.B_A2_ref = B_A2_ref
            self.B_B2_ref = B_B2_ref
            '''

        else:
            #st.write('h2cl2')
            m1 = self.massaH
            m2 = self.massaH
            m3 = self.massaCl
            m4 = self.massaCl
            B_ref = [-66.0103, -52.9951, -42.88, -34.849, -28.3552, -23.0219,  -18.5833,  -14.8471, -11.6713, -8.94904,  -6.59832, -4.55533, -2.76969, -1.20109,  0.18304, 1.40931,
                        2.49968,  3.47237,  4.34266, 5.12343, 5.82562, 6.45854,  7.0302, 7.54749, 8.01638,  8.44207, 8.82907, 9.18137, 9.50244,  9.79537,  10.0629,  10.3074,  10.531]
            T_ref = [300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700,
                        725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000, 1025,  1050,  1075, 1100 ]
            Req1 = 0.744013
            Req2 = 2.007880

            self.mi =  (m1 + m2) * (m3 + m4) / (m1 + m2 + m3 + m4) 
            #mi1 = (m1 + m2) / 2 * Req1**2 #/ 6.02252
            #mi2 = (m3 + m4) / 2 * Req2**2 #/ 6.02252

            self.r1 = (m2 * Req1) / (m1+m2) #(m2 * self.media_Reqs) / (m1 - m2)
            self.r2 = Req1 - self.r1
            self.r3 = (m4 * Req2) / (m3+m4) #(m4 * self.media_Reqs) / (m3 - m4)
            self.r4 = Req2 - self.r3

            self.I1 = m1 * self.r1**2 + m2 * self.r2**2
            self.I2 =  m3 * self.r3**2 + m4 * self.r4**2


            self.mi1 = (m1 * m2)/(m1 + m2)
            self.mi2 = (m3 * m4)/(m3 + m4)

            self.Ma = m1 + m2  # verify
            self.Mb = m3 + m4  # esses dados sao do h2cl2
            self.w1 = -0.216
            self.w2 = 0.279
            self.rho_1 = 0.026667
            self.rho_2 = 0.05087
            self.zc1 = 0.303
            self.zc2 = 0.073
            self.x1 = 0.5
            self.x2 = 0.5

            def B_A2_ref(T):
                return self.BCl2_ref(T)

            def B_B2_ref(T):
                return self.BH2_ref(T)
            self.B_A2_ref = B_A2_ref
            self.B_B2_ref = B_B2_ref
        #st.write('self.media_Reqs = ', self.media_Reqs, 'r1 = ',r1, 'r2 = ',r2, 'r3 = ', r3, 'r4 = ',r4 ,  'self.I1 = ', self.I1,  'self.I2 = ', self.I2,  'mi1 = ', mi1,  'mi2 = ', mi2)
            


    def calculate_virial_state(self):
        self.gen_ref()
        for T in range(60, 500, 5):
            B_cruzado = (0.25*(22.1*self.Mb + self.Ma*(22.1 + self.Mb*T)) * (0.083 + 0.0695*self.w1 + 0.0695*self.w2) * (
                self.rho_1**(1/3) + self.rho_2**(1/3))**3) / ((10.9*self.Ma + 10.9*self.Mb + self.Ma*self.Mb*T)*(self.zc1 + self.zc2) * self.rho_1 * self.rho_2)

            B_virial_state = 2*self.x1*self.x2*B_cruzado + \
                self.x2**2*self.B_A2_ref(T) + self.x1**2*self.B_B2_ref(T)

            self.B_virial_state_ref.append(B_virial_state)
            self.T_state.append(T)
        return self.B_virial_state_ref, self.T_state




class TemperatureSetter:
    def set_temperature(self, T):
        #st.write('t good')
        self.T = T

    

class Integrand:
    def __init__(self, um_final_instance, temperature_instance, ref_making_instance, derivatives_instance):
        self.derivatives_instance = derivatives_instance
        self.ref_making_instance = ref_making_instance
        self.ref_making_instance.gen_ref()
        
        self.um_final_instance = um_final_instance
        self.temperature_instance = temperature_instance
        self.mi = self.ref_making_instance.mi
        self.I1 = self.ref_making_instance.I1
        self.I2 = self.ref_making_instance.I2
        self.mi1 = self.ref_making_instance.mi1
        self.mi2 = self.ref_making_instance.mi2
        #st.write('integrand criada')

    def integrand_vegas(self, x):
        self.T = self.temperature_instance.T
        r, th_a, th_b, phi = x
        # Pass LC as an argument to UM_FINAL
        F = self.um_final_instance.calculate(r, th_a, th_b, phi) * rk / self.T
        if F < -1:
            F = -1

        #ff = constants_fixed.N_A / 4 * np.sin(th_a) * np.sin(th_b) * (r ** 2) * (1 - np.exp(-F))
        ff = 6.02252e-1 / 4 * np.sin(th_a) * np.sin(th_b) * (r ** 2) * (1 - np.exp(-F))
        return ff
    def integrand_c1(self, x):  # primeira correção quantica
        self.T = self.temperature_instance.T
        
        r, th_a, th_b, phi = x  
        F = self.um_final_instance.calculate(r, th_a, th_b, phi) * 7.24356 / self.T
        if F < -1:
            F = -1
        d_1r = self.derivatives_instance.d1r(r, th_a, th_b, phi)



        #c1 =  (10000* constants_fixed.N_A * constants_fixed.h**2 * constants_fixed.kb**3)/(48*self.mi*constants_fixed.uma*self.T**3) * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * d_1r**2 * r**2
        c1 = np.pi/12 * 0.602252 * 1.05459**2/(self.mi) * np.exp(-F) * 7.24356**3 / self.T**3 #* np.sin(th_a) * np.sin(th_b) * r**2 * d_1r**2
         
        return c1

    def integrand_c2(self,x):  # segunda correção quantica
        
        self.T = self.temperature_instance.T

        r, th_a, th_b, phi = x  
        F = self.um_final_instance.calculate(r, th_a, th_b, phi) * rk / self.T
        if F < -1:
            F = -1

        d_1th_a = self.derivatives_instance.d1th_a(r, th_a, th_b, phi)
        d_1th_b = self.derivatives_instance.d1th_b(r, th_a, th_b, phi)
        d_2th_a = self.derivatives_instance.d2th_a(r, th_a, th_b, phi)
        d_2th_b = self.derivatives_instance.d2th_b(r, th_a, th_b, phi)
        d_2phi = self.derivatives_instance.d2phi(r, th_a, th_b, phi)
        
        

        #VL1 = (constants_fixed.h**2)/(2*self.I1) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/np.sin(th_a)**2 * d_2phi)
        #VL2 = (constants_fixed.h**2)/(2*self.I2) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/np.sin(th_b)**2 * d_2phi)

        #c2 = -(constants_fixed.N_A * constants_fixed.kb**3)/(48*constants_fixed.uma*self.T**3) * np.exp(-F) * (VL1 + VL2) * r**2
        c2 = -np.pi/48 * 0.602252 * 1.05459**4/(self.mi) * np.exp(-F) * rk**4 / self.T**4
        return c2
    

    def integrand_c3(self, x):  # segunda correção quantica

        self.T = self.temperature_instance.T
        r, th_a, th_b, phi = x  
        F = self.um_final_instance.calculate(r, th_a, th_b, phi) * rk / self.T
        if F < -1:
            F = -1

        d_1th_a = self.derivatives_instance.d1th_a(r, th_a, th_b, phi)
        d_1th_b = self.derivatives_instance.d1th_b(r, th_a, th_b, phi)
        d_2th_a = self.derivatives_instance.d2th_a(r, th_a, th_b, phi)
        d_2th_b = self.derivatives_instance.d2th_b(r, th_a, th_b, phi)
        d_2phi = self.derivatives_instance.d2phi(r, th_a, th_b, phi)

        VL1 = (constants_fixed.h**2)/(2*self.mi1*r**2) * (1/np.tan(th_a) * d_1th_a + d_2th_a + 1/np.sin(th_a)**2 * d_2phi)
        VL2 = (constants_fixed.h**2)/(2*self.mi2*r**2) * (1/np.tan(th_b) * d_1th_b + d_2th_b + 1/np.sin(th_b)**2 * d_2phi)


        #c3 = -(constants_fixed.N_A * constants_fixed.kb**3)/(48*constants_fixed.uma*self.T**3) * np.exp(-F) * (VL1 + VL2) * r**2
        c3_a = -np.pi/24 * 0.602252 * 1.05459**2/(self.I1) * np.exp(-F) * rk**2 / self.T**2
        c3_b = -np.pi/24 * 0.602252 * 1.05459**2/(self.I2) * np.exp(-F) * rk**2 / self.T**2
        return c3_a + c3_b
    

    def integrand_c4(self, x):  # quarta correção quantica

        self.T = self.temperature_instance.T
        r, th_a, th_b, phi = x  
        F = self.um_final_instance.calculate(r, th_a, th_b, phi) * rk / self.T
        if F < -1:
            F = -1

        d_1r = self.derivatives_instance.d1r(r, th_a, th_b, phi)
        d_2r = self.derivatives_instance.d2r(r, th_a, th_b, phi)
        d_1th_a = self.derivatives_instance.d1th_a(r, th_a, th_b, phi)
        d_1th_b = self.derivatives_instance.d1th_b(r, th_a, th_b, phi)
        d_2th_a = self.derivatives_instance.d2th_a(r, th_a, th_b, phi)
        d_2th_b = self.derivatives_instance.d2th_b(r, th_a, th_b, phi)
        d_2phi = self.derivatives_instance.d2phi(r, th_a, th_b, phi)



        #c4 = -(constants_fixed.N_A * constants_fixed.h**4 * constants_fixed.kb**4)/(1920 * constants_fixed.uma*self.mi**2 * self.T**4) * np.sin(th_a) * np.sin(th_b) * np.exp(-F) * \
        #(d_2r**2  + 2/r**2 * d_1r**2 + (10 * constants_fixed.kb)/(9*self.T*r) * d_1r**3 - (5* constants_fixed.kb**2)/(36*self.T**2) * d_1r**4) * r**2                                      

        c4 = -np.pi/24 * 0.602252 * 1.05459**2/(self.mi) * np.exp(-F) * rk**2 / self.T**2

        #c4a = -N_A * h**4 * rk**4 / (1920 * mi**2 * T**4) * np.sin(th_a) * np.sin(th_b) * np.exp(-F)
        #c4b = (d_2r**2 + 2/r**2*(d_1r**2) + 10 * rk/(9 * T * r) * d_1r**3 - 5 * rk**2/(36 * T**2) * d_1r**4) * r**2
        #c4aa = c4a * c4b
        # c4p = np.sin(th_a) * np.sin(th_b) * np.exp(-F) * (d_2r**2 + 2/r**2*(d_1r**2) + 10 * rk/(9 * T * r) * d_1r**3 - 5 * rk**2/(36 * T**2) * d_1r**4) * r**2  #talvez sem rk/T
        #c4 = -np.pi/12 * N_A * h**2/mi * c4p *(rk/T)**2
        return c4
    





class Mean_req:
    def __init__(self, data):
        self.data = data
        #st.write('mean_req criada')

    def calculate_mean(self):
        Reqs = [self.data.loc[0, 'Req'], self.data.loc[1, 'Req'], self.data.loc[2, 'Req'], 
                self.data.loc[3, 'Req'], self.data.loc[0, 'Req'], self.data.loc[0, 'Req']]

        return mean(Reqs)


class Calculate_virialis:
    def __init__(self, integrand_instance, mean_req_instance, temperature_instance):
        self.integrand_instance = integrand_instance
        self.mean_req_instance = mean_req_instance
        self.temperature_instance = temperature_instance
        #st.write('calculate_virial criada')
        self.B_clas = []
        self.Tstr = []
        self.B_main = []


    def compute_results_main(self):

        #st.write('outisde look good')
        # Integrator and definition of integration limits
        integ = vegas.Integrator(
        [[0.2 * self.mean_req_instance.calculate_mean(), 2 * self.mean_req_instance.calculate_mean()],
        [0, np.pi], [0, np.pi], [0, 2 * np.pi]])

        #logging.info(f"Starting integration for temperature {temp}")
        #st.write(f"Starting integration for temperature {temp}")

        self.temperature_instance.set_temperature(temp)
        result = integ(self.integrand_instance.integrand_vegas, nitn=10, neval=10000)
        st.write('result class = ', result)
        result_c1  = integ(self.integrand_instance.integrand_c1, nitn=10, neval=10000)
        st.write('result_c1 = ', result_c1)
        result_c2  = integ(self.integrand_instance.integrand_c2, nitn=10, neval=10000)
        st.write('result_c2 = ', result_c2)
        result_c3  = integ(self.integrand_instance.integrand_c3, nitn=10, neval=10000)
        st.write('result_c3 = ', result_c3)
        result_c4  = integ(self.integrand_instance.integrand_c4, nitn=10, neval=10000)
        st.write('result_c4 = ', result_c4)
        
        self.B_clas.append(result.mean)
        self.Tstr.append(temp)
        self.B_main.append(result.mean + result_c1.mean + result_c2.mean + result_c3.mean + result_c4.mean)

        return self.B_clas, self.Tstr, self.B_main
        
    
    
class graph_gen:

    def __init__(self, LCs_instance, momentum_instance, B_values, T_values, B_main_values, T_state_mix, B_state_mix, LCs_instance_ryd = None, momentum_instance_ryd = None, B_values_ryd = None, T_values_ryd = None, B_main_values_ryd = None):
        self.B_values = B_values
        self.T_values = T_values
        self.B_main_values = B_main_values
        self.T_state_mix = T_state_mix
        self.B_state_mix = B_state_mix
        self.r = np.linspace(0, 10, 100)
        self.th_a = np.linspace(0, np.pi, 100)
        self.th_b = np.linspace(0, np.pi, 100)
        self.phi = np.linspace(0, 2*np.pi, 100)
        #st.write('graph_gen criada')
        self.LCs_instance = LCs_instance
        self.momentum_instance = momentum_instance

        self.LCs_instance_ryd = LCs_instance_ryd
        self.momentum_instance_ryd = momentum_instance_ryd

        self.B_values_ryd = B_values_ryd
        self.T_values_ryd = T_values_ryd
        self.B_main_values_ryd = B_main_values_ryd



    def plot_LCs(self):

        self.UH = self.LCs_instance.UH(self.r)
        self.UL = self.LCs_instance.UL(self.r)
        self.UX = self.LCs_instance.UX(self.r)
        self.UZ = self.LCs_instance.UZ(self.r)
        self.UTa = self.LCs_instance.UTa(self.r)
        self.UTb = self.LCs_instance.UTb(self.r)
        self.De_lim = self.LCs_instance.imported.h_De
        
        fig_LCs = go.Figure()

        fig_LCs.add_trace(
            go.Scatter(x=self.r, y=self.UH, name='UH')
        )
        fig_LCs.add_trace(
            go.Scatter(x=self.r, y=self.UX, name='UX')
        )
        fig_LCs.add_trace(
            go.Scatter(x=self.r, y=self.UZ, name='UZ')
        )
        fig_LCs.add_trace(
            go.Scatter(x=self.r, y=self.UTa, name='UTa')
        )
        fig_LCs.add_trace(
            go.Scatter(x=self.r, y=self.UTb, name='UTb')
        )
        fig_LCs.add_trace(
            go.Scatter(x=self.r, y=self.UL, name='UL')
        )

        fig_LCs.update_xaxes(range=[0, 10])
        fig_LCs.update_yaxes(range=[-3*self.De_lim, 3*self.De_lim])


        fig_LCs.update_layout(height=600, width=800,
                        title_text="Graph of Leading Configurations per Distance", xaxis_title=r"$r[A]$",
                        yaxis_title='U[kcal/mol]'
                            )
        st.plotly_chart(fig_LCs, use_container_width=True)

    def plot_harm_esf(self):

        self.UM_000 = self.momentum_instance.UM_000(self.r)
        self.UM_202 = self.momentum_instance.UM_202(self.r)
        self.UM_022 = self.momentum_instance.UM_022(self.r)
        self.UM_220 = self.momentum_instance.UM_220(self.r)
        self.UM_222 = self.momentum_instance.UM_222(self.r)
        self.UM_224 = self.momentum_instance.UM_224(self.r)
        self.De_lim = self.LCs_instance.imported.h_De
        
        fig_harm_esf = go.Figure()
        fig_harm_esf.add_trace(
            go.Scatter(x=self.r, y=self.UM_000, name='UM_000')
        )
        fig_harm_esf.add_trace(
            go.Scatter(x=self.r, y=self.UM_202, name='UM_202')
        )
        fig_harm_esf.add_trace(
            go.Scatter(x=self.r, y=self.UM_220, name='UM_220')
        )
        fig_harm_esf.add_trace(
            go.Scatter(x=self.r, y=self.UM_022, name='UM_022')
        )
        fig_harm_esf.add_trace(
            go.Scatter(x=self.r, y=self.UM_222, name='UM_222')
        )
        fig_harm_esf.add_trace(
            go.Scatter(x=self.r, y=self.UM_224, name='UM_224')
        )
        fig_harm_esf.update_xaxes(range=[0, 10])
        fig_harm_esf.update_yaxes(range=[-3*self.De_lim, 5*self.De_lim])
        fig_harm_esf.update_layout(height=600, width=800,
                        title_text="Graph of Complex Energies per Distance", xaxis_title=r"$r[A]$",
                        yaxis_title='U[kcal/mol]')
        st.plotly_chart(fig_harm_esf, use_container_width=True)

    def plot_virial(self):
        fig_virial = go.Figure()


        fig_virial.add_trace(
            go.Scatter(x=self.T_values, y=self.B_values, name='B classic')
        )

        fig_virial.add_trace(
            go.Scatter(x=self.T_values, y=self.B_main_values, name='B with corrections')
        )

        fig_virial.add_trace(
            go.Scatter(x=self.T_state_mix, y=self.B_state_mix,
                        name='B virial state reference')
        )


        fig_virial.update_layout(height=600, width=800,
        title_text="Graph of the Second Virial Coefficient B per Temperature T", yaxis_title='B(T)[cm^3/mol]',
        xaxis_title="T[K]")

        st.plotly_chart(fig_virial, use_container_width=True)
    
    def plot_virial_both(self):
        fig_virial = go.Figure()

        fig_virial.add_trace(
            go.Scatter(x=self.T_values_ryd, y=self.B_main_values_ryd, name='B rydberg')
        )

        fig_virial.add_trace(
            go.Scatter(x=self.T_values, y=self.B_main_values, name='B ILJ')
        )

        fig_virial.add_trace(
            go.Scatter(x=self.T_state_mix, y=self.B_state_mix,
                        name='B virial state reference')
        )


        fig_virial.update_layout(height=600, width=800,
        title_text="Graph of the Second Virial Coefficient B per Temperature T", yaxis_title='B(T)[cm^3/mol]',
        xaxis_title="T[K]")

        st.plotly_chart(fig_virial, use_container_width=True)
            
if st.button('Calculate'):
    if uploaded_file is not None:
        
        if potential == 'Improved Leonard-Jonnes Potential':
            imported = Data_manage(uploaded_file, st)
            imported.process_data_ilj()
            if imported.data is not None:  # Check if data has been processed
                LCs_instance = LCs_ilj(imported)
                momentum_instance = momentum_ilj()
                derivatives_instance = derivatives_ilj()

                #st.write('momentum', stop_momentum - start)

                um_final_instance = UM_FINAL(momentum_instance)
                #st.write('um_final', stop_umfinal - start)

                temperature_instance = TemperatureSetter()
                mean_req_instance = Mean_req(imported.data)
                media_Reqs = mean_req_instance.calculate_mean()
                ref_making_instance = RefMaking(gas, media_Reqs, temperature_instance)
                B_state_mix, T_state_mix = ref_making_instance.calculate_virial_state()
                integrand_instance = Integrand(um_final_instance, temperature_instance, ref_making_instance, derivatives_instance)
                #st.write('integrand', stop_integrand - start)
                

                
                virial_calculation_instance = Calculate_virialis(integrand_instance, mean_req_instance,temperature_instance)
                #st.write('virial_calc', stop_virial - start)

                for temp in range(50, 1000, step):
                    st.write('T=', temp)

                    B_values, T_values, B_main_values = virial_calculation_instance.compute_results_main()
                #st.write('LC', stop_calc_instance - start)


                #st.write("B values:", B_values)
                #st.write("T values:", T_values)
                #st.write("Main B values:", B_main_values)
                #st.write("Main B values:", B_main_values)


                graph_instance = graph_gen(LCs_instance, momentum_instance, B_values, T_values, B_main_values, T_state_mix, B_state_mix)
                graph_instance.plot_LCs()
                graph_instance.plot_harm_esf()
                graph_instance.plot_virial()



        elif potential == 'Rydberg Potential':
            imported = Data_manage(uploaded_file, st)
            imported.process_data_ryd()
            if imported.data is not None:  # Check if data has been processed
                LCs_instance = LCs_ryd(imported)
                momentum_instance = momentum_ryd()
                derivatives_instance = derivatives_ryd()

                #st.write('momentum', stop_momentum - start)

                um_final_instance = UM_FINAL(momentum_instance)
                #st.write('um_final', stop_umfinal - start)

                temperature_instance = TemperatureSetter()
                mean_req_instance = Mean_req(imported.data)
                media_Reqs = mean_req_instance.calculate_mean()
                ref_making_instance = RefMaking(gas, media_Reqs, temperature_instance)
                B_state_mix, T_state_mix = ref_making_instance.calculate_virial_state()
                integrand_instance = Integrand(um_final_instance, temperature_instance, ref_making_instance, derivatives_instance)
                #st.write('integrand', stop_integrand - start)
                

                
                virial_calculation_instance = Calculate_virialis(integrand_instance, mean_req_instance,temperature_instance)
                #st.write('virial_calc', stop_virial - start)

                for temp in range(50, 1000, step):
                    st.write('T=', temp)

                    B_values, T_values, B_main_values = virial_calculation_instance.compute_results_main()
                #st.write('LC', stop_calc_instance - start)


                #st.write("B values:", B_values)
                #st.write("T values:", T_values)
                #st.write("Main B values:", B_main_values)
                #st.write("Main B values:", B_main_values)


                graph_instance = graph_gen(LCs_instance, momentum_instance, B_values, T_values, B_main_values, T_state_mix, B_state_mix)
                graph_instance.plot_LCs()
                graph_instance.plot_harm_esf()
                graph_instance.plot_virial()
        elif potential == "both":
            imported = Data_manage(uploaded_file, st)
            imported.process_data_ryd()
            imported.process_data_ilj()
            
            if imported.data_ilj is not None and imported.data_ryd is not None:  # Check if data has been processed

                LCs_instance_ryd = LCs_ryd(imported)
                LCs_instance_ilj = LCs_ilj(imported)

                momentum_instance_ryd = momentum_ryd()
                momentum_instance_ilj = momentum_ilj()

                derivatives_instance_ryd = derivatives_ryd()
                derivatives_instance_ilj = derivatives_ilj()


                #st.write('momentum', stop_momentum - start)

                um_final_instance_ryd = UM_FINAL(momentum_instance_ryd)
                um_final_instance_ilj = UM_FINAL(momentum_instance_ilj)
                #st.write('um_final', stop_umfinal - start)

                temperature_instance = TemperatureSetter()
                mean_req_instance_ilj = Mean_req(imported.data_ilj)
                mean_req_instance_ryd = Mean_req(imported.data_ryd)
                media_Reqs_ilj = mean_req_instance_ilj.calculate_mean()
                media_Reqs_ryd = mean_req_instance_ryd.calculate_mean()
                ref_making_instance = RefMaking(gas, media_Reqs_ilj, temperature_instance)

                B_state_mix, T_state_mix = ref_making_instance.calculate_virial_state()
                integrand_instance_ryd = Integrand(um_final_instance_ryd, temperature_instance, ref_making_instance, derivatives_instance_ryd)
                integrand_instance_ilj = Integrand(um_final_instance_ilj, temperature_instance, ref_making_instance, derivatives_instance_ilj)
                #st.write('integrand', stop_integrand - start)
                
                virial_calculation_instance_ryd = Calculate_virialis(integrand_instance_ryd, mean_req_instance_ryd,temperature_instance)

                for temp in range(50, 1000, step):
                    st.write('T=', temp)
                    st.write('ryd calculating...')
                    B_values_ryd, T_values_ryd, B_main_values_ryd = virial_calculation_instance_ryd.compute_results_main()
                    st.write(B_main_values_ryd, 'ryd')
                #st.write('LC', stop_calc_instance - start)
                del virial_calculation_instance_ryd

                virial_calculation_instance_ilj = Calculate_virialis(integrand_instance_ilj, mean_req_instance_ilj,temperature_instance)
                #st.write('virial_calc', stop_virial - start)

                for temp in range(50, 1000, step):
                    st.write('T=', temp)
                    st.write('ILJ calculating...')
                    B_values_ilj, T_values_ilj, B_main_values_ilj = virial_calculation_instance_ilj.compute_results_main()
                    st.write(B_main_values_ilj, 'ilj')
                #st.write('LC', stop_calc_instance - start)

                del virial_calculation_instance_ilj
                #st.write("B values:", B_values)
                #st.write("T values:", T_values)
                #st.write("Main B values:", B_main_values)
                #st.write("Main B values:", B_main_values)

                graph_instance = graph_gen(LCs_instance_ilj, momentum_instance_ilj, B_values_ilj, T_values_ilj, B_main_values_ilj, T_state_mix, B_state_mix, LCs_instance_ryd, momentum_instance_ryd, B_values_ryd, T_values_ryd, B_main_values_ryd)
                graph_instance.plot_LCs()
                graph_instance.plot_harm_esf()
                graph_instance.plot_virial_both()
            


       
    else:
        st.info('☝️ Upload a .dat file')
