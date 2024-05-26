# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 18:27:59 2021

@author: Alberto
"""



#IMPORTAÇÃO DAS BIBLIOTECAS

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
from scipy.integrate import quad
import random
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go


# TIME MEASURE

start = time.time()



#ATRIBUIÇÃO DOS DADOS DE ENTRADA A CONSTANTES
#atribuindo os coeficientes para cada uma das linhas do arquivo aberto



st.title('Virial Coefficient of a Diatomic Molecule AB')

uploaded_file = st.file_uploader("upload a file")


potential = st.selectbox(
    'What potential energy do you want to use?',
    ('Rydberg Potential', 'Improved Leonard-Jonnes Potential'))



step = st.selectbox(
    'What is the Temperature step that you want to use?',
    (100, 50, 25, 200))





st.write('The calculus will be made using ', potential,
         'on a gas of structure AB, using a temperature step of', step)



atom1 = st.selectbox(
    'What is your first atom?',
    ('H', 'F', 'Cl', 'Br', 'O', 'C'))
atom2 = st.selectbox(
    'What is your second atom?',
    ('H', 'F', 'Cl', 'Br','O'))
gas = atom1 + atom2



#DEFININDO CONSTANTES
kb = 0.0019872041*1000# para kcal/mol, *1000 nao sei pq //#8.617332478e-5 eV/K
rk= 1/(kb)    #4.184/(1.38064853*6.02252) #   503.188016899  #7.24297 #*10**22
h=1.05451  #*10**-34
N_A=0.602252
r0=  1.33 #Req
irk= 1.0/rk
K=-272.15

#DERIVADAS NUMÉRICAS QUE SERÃO USADAS

def derivative(g,a,method='central',p=0.01): #derivada de promeira ordem
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
        
def monte_carlo_uniform(func, a=0, b=1, n=1000):
    subsets = np.arange(0,n+1,n/10)
    steps=n/10
    u=np.zeros(n)
    for i in range(10):
        start = int(subsets[i])
        end = int(subsets[i+1])
        u[start:end] = np.random.uniform(low=i/10, high=(i+1)/10, size=end-start)
    np.random.shuffle(u)
    u_func=func(a+(b-a)*u)
    s= ((b-a)/n)*u_func.sum()
    
    return s            
     

if st.button('Calculate'):
    if uploaded_file is not None:

        if potential == 'Improved Leonard-Jonnes Potential' :
            st.write('LENNARD-JONNES IS IN CONSTRUCTION')
            data = pd.read_csv(uploaded_file, sep="\s+", header=None)
            data.columns = ["alpha", "beta", "mp", "De", "Req" ]
            st.subheader('DataFrame')
            st.write(data)


            alpha = float(data['alpha'])
            beta = float(data['beta'])
            mp = float(data['mp'])
            De = float(data['De'])
            Req = float(data['Req'])
            
            def U(r):

                n = beta + alpha * (r / Req) ** 2
                return De * ((mp/(n - mp) * (Req/r) ** n) - (n/(n - mp) * (Req/r) ** mp))

            def fprincipal(r):
                F=U(r)*rk/T

                return (r**2)*(1-np.exp(-F))#*10**-8

        elif potential == 'Rydberg Potential':
        
            data = pd.read_csv(uploaded_file, sep="\s+", header=None)
            data.columns = ["a1", "a2", "a3", "a4", 'a5', 'De', 'Req', 'Eref']
            st.subheader('DataFrame')
            st.write(data)

            a1 = float(data['a1'])
            a2 = float(data['a2'])
            a3 = float(data['a3'])
            a4 = float(data['a4'])
            a5 = float(data['a5'])
            De = float(data['De'])
            Req = float(data['Req'])
            Eref = float(data['Eref'])



            print(a1, a2, De, Req)


            #POTENCIAL
            def U(r):
                y = (r-Req)/Req
                result = -De*(1 + a1*y + a2*y**2 + a3*y**3 +
                a4*y**4 + a5*y**5) * np.exp(-a1*y) + Eref
                return result
                

            #FUNÇÃO PRINCIPAL QUE SERÁ INTEGRADA
            def fprincipal(r):
                F=U(r)*rk/T

            
                return (r**2)*(1-np.exp(-F))

        




#INPUT DAS MASSAS PARA O MI


        B_ref = []
        if atom1=='O':
            atom1=15.99491
            if atom2 == 'O':
                T_ref = [70  ,85 ,100   ,120 ,140,170 ,200,240 ,280,330 ,380  ,430,495]
                B_ref =[ -369.9,-262.9 ,-195.2,-138.8  ,-103.5  ,-69.7, -49.9, -32.2, -20.3, -10.1 , -2.9 , 2.6 ,8.0] 
            
        elif atom1=='F':
            atom1= 18.99840
            if atom2 == 'F':
                T_ref = [85,95,105,120,135,150,170,190,210,230,250,270,295]
                B_ref = [-212.6, -172.1 ,-142.4 ,-110.3 ,-87.7 ,-71.0 ,-54.6 ,-42.5 ,-33.2 ,-25.9 ,-20.1 ,-15.2 ,-10.3]
                
        elif atom1=='C':
            atom1= 12.0
            if atom2 == 'O':
                T_ref = [125,135,145,165,185,205,235,265,295,325,355,395,435,475,515,555,570]
                B_ref = [-118.1,-101.3,-87.6,-66.5,-51.2,-39.4,-26.3,-16.6,-9.1,-3.2,1.6,6.8,10.9,14.3,17.1,19.5,20.3]
        elif atom1=='N':
            atom1= 14.00307
            
            if atom2 == 'O':
                T_ref = [125,135,150,165,185,215,245,285,325,375,425,475]
                B_ref = [ -201.5, -158.1, -115.9, -89.4  , -67.1, -47.5, -35.2, -23.8, -15.1, -6.2, 1.3, 7.8]
            elif atom2 == 'N':
                T_ref = [75,85,100,120,140,170,200,240,280,330,380,440,500,570,640,745]
                B_ref = [ -276.8, -218.0, -160.7, -113.6 , -83.4, -54.4, -35.9, -19.6, -8.8, 0.5, 6.9, 12.4, 16.4, 19.8,22.5,25.3]
        elif atom1=='H':
            atom1= 1.00783
            if atom2 == 'H':
                T2_ref = [15,20,25,30,35,40,45,50,55,60,65,75,85,100,115,130,145,165,180,200,220,250,290,330,370,410,450,490]
                T_ref= x = np.arange(60, 500, 10)
                B_ref= 1.7472*10 - (1.2926* 10**2)/(x) - (2.6988* 10**5)/(x)**2 + (8.0282 * 10**6)/(x)**3

            


            B2_ref = [-239.2,-150.6,-105.7,-79.0,-61.4,-49.0,-39.8,-32.7,-27.1,-22.6,-19.2,-13.2,-8.3,-2.8,1.2,4.2,6.4,8.6,9.8,11.1,12.1,13.2,14.1,14.8,15.3,15.7,15.9,16.2]

        if atom2=='O':
            atom2=15.99491
        elif atom2=='F':
            atom2= 18.99840
        elif atom2=='C':
            atom2= 12.0
        elif atom2=='N':
            atom2= 14.00307
        elif atom2=='H':
            atom2= 1.00783
            
            
        '''  
        massaC = 12.0
        massaO=15.99491
        massaN=14.00307
        massaH=1.00783
        massaF=18.99840 
        massaBr = 79.904
        massaCl = 34.96885
        m1=massaC
        m2=massaO
        '''

        mi = ((atom1*atom2))/(atom1+atom2) #definição do mi 0.602252

        #cria as listas que serão preenchidas para fazer o grafico
        Bstr=[]
        Tstr=[]
        Bmonte=[]

        R_0 = 1.15
        #CALCULO DAS INTEGRAIS E COEFICIENTES VIRIAIS PARA DIFERENTES TEMPERATURAS
        for T in range(10, 100, 20): #CALCULOS QUE SERÃO FEITOS NO PRIMEIRO INTERVALO DE TEMPERATURA
        
            integralprinc1, err0i1 = quad(fprincipal,0,np.inf) #integração


            
            B1=(2*np.pi*N_A)*(integralprinc1+(R_0**3)/3)# #coef do virial principal 


            inte1=monte_carlo_uniform(fprincipal,a=0,b=10,n=1800)
            B1mc=(2*np.pi*N_A)*(inte1+(R_0**3)/3)    
            
            
        

            Btotal1= (B1) #soma de cada um dos coef encontrados resultando no coeficiente total
            


            
            Tstr.append(T) #temperaturas usadas são colocadas na lista
            
            Bstr.append(Btotal1) #Coef total achado é alocado na lista

            Bmonte.append(B1mc)
            print('B=', Btotal1, 'para T=',T)
            
        for T in range(100, 300, 35):  #CALCULOS QUE SERÃO FEITOS NO SEGUNDO INTERVALO DE TEMPERATURA

            integralprinc2, err0i2 =quad(fprincipal,0,np.inf) #integração
        # print(integral, err3)
            B2=(2*np.pi*N_A)*(integralprinc2+(R_0**3)/3)
        
         
            
            Btotal2= (B2) #soma de cada um dos coef encontrados resultando no coeficiente total
            inte2=monte_carlo_uniform(fprincipal,a=0,b=10,n=1800)
            B2mc=(2*np.pi*N_A)*(inte1+(R_0**3)/3)              
            #d1 = decimal.Decimal(B1)
            
            Tstr.append(T) #temperaturas usadas são colocadas na lista
            

            Bstr.append(Btotal2) #Coef total achado é alocado na lista
            Bmonte.append(B2mc)

            print('B=',Btotal2, 'para T=',T)
            
            
            
        for T in range(300, 500, 70):  #CALCULOS QUE SERÃO FEITOS NO TERCEIRO INTERVALO DE TEMPERATURA
        


            integralprinc3, err0i3 =quad(fprincipal,0,np.inf) #integração
            #print(integral, err3)
            B3=(2*np.pi*N_A)*(integralprinc3+(R_0**3)/3) #coef do virial principal
            
            inte3=monte_carlo_uniform(fprincipal,a=0,b=10,n=1800)
            B3mc=(2*np.pi*N_A)*(inte1+(R_0**3)/3)    
            
            Btotal3= (B3)  #soma de cada um dos coef encontrados resultando no coeficiente total
            
            #d1 = decimal.Decimal(B1)
            
            Tstr.append(T) #temperaturas usadas são colocadas na lista
            
            Bstr.append(Btotal3) #Coef total achado é alocado na lista
            Bmonte.append(B3mc)

            print('B=',Btotal3, 'para T=',T)



        r = np.linspace(0, 10, 100)

     
        Ur = U(r)
      
            
   
        st.write(Ur, type(Ur))



        #CRIANDO O GRÁFICO


        # Create two subplots and unpack the output array immediately

        fig1 = go.Figure()

        fig1.add_trace(
            go.Scatter(x=r, y=Ur, name='SEP'))
        
        fig1.update_layout(height=600, width=800,
                            title_text="Graph of the PES", yaxis_title='U(eV)',
                            xaxis_title="R")




        fig2 = go.Figure()

        fig2.add_trace(
            go.Scatter(x=Tstr, y=Bstr, name='dados calculados'))
        fig2.add_trace(
            go.Scatter(x=Tstr, y=B_ref, name='dados de referencia'))
        fig2.add_trace(
            go.Scatter(x=Tstr, y=Bmonte, name='dados calculados com Monte Carlo'))


        fig2.update_layout(height=600, width=800,
                            title_text="Graph of the Second Virial Coefficient B per Temperature T", yaxis_title='B(T)[cm^3/mol]',
                            xaxis_title="T[K]")
        #fig2.update_yaxes(range=[-1000000, 5000000])


        st.plotly_chart(fig1, use_container_width=True)

        st.plotly_chart(fig2, use_container_width=True)



        st.success('Calculations finished!', icon="✅")