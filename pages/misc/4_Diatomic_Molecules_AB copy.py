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
    ('H', 'F', 'Cl', 'Br', 'O'))
atom2 = st.selectbox(
    'What is your second atom?',
    ('H', 'F', 'Cl', 'Br','O'))
gas = atom1 + atom2



#DEFININDO CONSTANTES

rk=4.184/(1.38064853*6.02252) #   503.188016899  #7.24297 #*10**22
h=1.05451  #*10**-34
N_A=6.02252
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
        
        
     

if st.button('Calculate'):
    if uploaded_file is not None:

        if potential == 'Improved Leonard-Jonnes Potential' :
                st.write('LEONARD-JONNES IS IN CONSTRUCTION')

        elif potential == 'Rydberg Potential':
        
            data = pd.read_csv(uploaded_file, sep="\s+", header=None)
            data.rows = ["a1", "a2", "a3", "a4", 'a5', 'De', 'Req', 'Eref']
            st.subheader('DataFrame')
            st.write(data)
            a1 = data.iloc[0]
            a2 = data.iloc[1]
            a3 = data.iloc[2]
            a4 = data.iloc[3]
            a5 = data.iloc[4]
            De = data.iloc[5]
            Req = data.iloc[6]
            Eref = data.iloc[7]


            print(a1, a2, De, Req)


            #POTENCIAL
            def U(r):
                y = (r-Req)/Req
                U = -De*(1 + a1*y + a2*y**2 + a3*y**3 +
                a4*y**4 + a5*y**5) * np.exp(-a1*y) + Eref
                return U
                

            #FUNÇÃO PRINCIPAL QUE SERÁ INTEGRADA
            def fprincipal(r):
                F=U(r)*rk/T

            
                return (r**2)*(1-np.exp(-F))#*10**-8

        




            def lim_inf(r):
                F=U(r)*rk/10

                
                return (r**2)*(1-np.exp(-F))#*10**-8


            def lim_sup(r):
                F=U(r)*rk/500

            
                return (r**2)*(1-np.exp(-F))#*10**-8




            def f1correção(r): #função usada na primeira correção
                
                    
                derivadinha=derivative(U,r,method='central',p=1e-8)
                
            
                        
                    #derivadinha=derivative(derivada, 1.0, dx=1e-6)
                    
                    #print(derivadinha)
                
                return (r**2)*(derivadinha**2)*(np.exp(-U(r)/(irk*T)))


            def f2correção(r): #função usada na segunda correção
                
            
                    
            
            
                        
                    #derivadinha=derivative(derivada, 1.0, dx=1e-6)
                    
                    #print(derivadinha)
                
                return (r**2)*(np.exp(-U(r)/(irk*T)))#*derivada parcial Dv/Dr

            def f3correção(r): #função usada na terceira correção
                
            
                    
            
            
                        
                    #derivadinha=derivative(derivada, 1.0, dx=1e-6)
                    
                    #print(derivadinha)
                
                return (r**2)*(np.exp(-U(r)/(irk*T)))#*derivada parcial Dv/Dr

            def f4correção(r): #função usada na quarta correção
                
            
                    
            
                
                    
                derivadinha=derivative(U,r,method='central',p=1e-8)
                derivadinha2=derivative2(U,r,method='central',p=1e-8)
                        
                    
                
                return (np.exp(-U(r)/(irk*T)))*(((derivadinha2)**2) + (2/(r**2))*(derivadinha**2) + (10/(9*irk*T*r))*((derivadinha)**3) - (5/(36*(irk**2)*(T**2)))*((derivadinha)**4)) * (r**2)#*derivada parcial Dv/Dr

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

        rstr=[]
        upot_str=[]
        #CALCULO DAS INTEGRAIS E COEFICIENTES VIRIAIS PARA DIFERENTES TEMPERATURAS
        for T in range(10, 100, 5): #CALCULOS QUE SERÃO FEITOS NO PRIMEIRO INTERVALO DE TEMPERATURA
        
            integralprinc1, err0i1 = quad(fprincipal,0,np.inf) #integração


            
            B1=(2*np.pi*N_A)*(integralprinc1+(r0**3)/3)# #coef do virial principal 

            
            Bmonte1=(2*np.pi*N_A)*(inte1+(r0**3)/3)
            
            corr1i1, err1i1 =quad(f1correção,0,np.inf) #integração
        
            B1corr1=((N_A*10**-24*(h**2))*(corr1i1))/(48*(irk**3)*(T**3)*mi) #coef do virial da primeira correção



            corr2i1, err2i1 =quad(f2correção,0,np.inf) #integração
            
            B1corr2=((-N_A*10**-24)*(corr2i1))/(48*(irk**3)*(T**3))  #coef do virial da segunda correção
            
        
            # d1 = decimal.Decimal(B1)
        
            corr3i1, err3i1 =quad(f3correção,0,np.inf) #integração
            #print(integral1)
            #print(integral1, T)
            B1corr3=((-N_A*10**-24)*(corr3i1))/(48*(irk**3)*(T**3))

        
            corr4i1, err4i1 =quad(f4correção,0,np.inf) #integração
        
        
            B1corr4=((-N_A*10**-24*(h**4))*(corr4i1))/(1920*(irk**4)*(T**4)*(mi**2)) #*mi
            
            
            Btotal1= (B1 + B1corr1 + B1corr2 + B1corr3 + B1corr4) #soma de cada um dos coef encontrados resultando no coeficiente total
            


            
            Tstr.append(T) #temperaturas usadas são colocadas na lista
            
            Bmonte.append(Bmonte1)
            Bstr.append(Btotal1) #Coef total achado é alocado na lista
            print('B=', Btotal1, 'para T=',T)
            
        for T in range(100, 300, 25):  #CALCULOS QUE SERÃO FEITOS NO SEGUNDO INTERVALO DE TEMPERATURA

            integralprinc2, err0i2 =quad(fprincipal,0,np.inf) #integração
        # print(integral, err3)
            B2=(2*np.pi*N_A)*(integralprinc2+(r0**3)/3)
        
            inte2=monte_carlo_uniform(fprincipal,a=0,b=2,n=1800)
            Bmonte2=(2*np.pi*N_A)*(inte2+(r0**3)/3)
            
            corr1i2, err1i2 =quad(f1correção,0,np.inf) #integração
        
            B2corr1=((N_A*10**-24*(h**2))*(corr1i2))/(48*(irk**3)*(T**3)*mi) #coef do virial da primeira correção



            corr2i2, err2i2 =quad(f2correção,0,np.inf) #integração
        
            B2corr2=((-N_A*10**-24)*(corr2i2))/(48*(irk**3)*(T**3)) #coef do virial da segunda correção
            
        

        
            corr3i2, err3i2 =quad(f3correção,0,np.inf) #integração

            B2corr3=((-N_A*10**-24)*(corr3i2))/(48*(irk**3)*(T**3)) #coef do virial da terceira correção

        
            corr4i2, err4i2 =quad(f4correção,0,np.inf) #integração
        
            B2corr4=((-N_A*10**-24*(h**4))*(corr4i2))/(1920*(irk**4)*(T**4)*(mi**2)) #coef do virial da quarta correção
            
            
            Btotal2= (B2 + B2corr1 + B2corr2 + B2corr3 + B2corr4) #soma de cada um dos coef encontrados resultando no coeficiente total
            
            #d1 = decimal.Decimal(B1)
            
            Tstr.append(T) #temperaturas usadas são colocadas na lista
            
            Bmonte.append(Bmonte2)

            Bstr.append(Btotal2) #Coef total achado é alocado na lista
            print('B=',Btotal2, 'para T=',T)
            
            
            
        for T in range(300, 500, 50):  #CALCULOS QUE SERÃO FEITOS NO TERCEIRO INTERVALO DE TEMPERATURA
        


            integralprinc3, err0i3 =quad(fprincipal,0,np.inf) #integração
            #print(integral, err3)
            B3=(2*np.pi*N_A)*(integralprinc3+(r0**3)/3) #coef do virial principal
            
            inte3=monte_carlo_uniform(fprincipal,a=0,b=2,n=1800)
            Bmonte3=(2*np.pi*N_A)*(inte3+(r0**3)/3)
        
            corr1i3, err1i3 =quad(f1correção,0,np.inf) #integração
        
            B3corr1=((N_A*10**-24*(h**2))*(corr1i3))/(48*(irk**3)*(T**3)*mi) #coef do virial da segunda correção



            corr2i3, err2i3 =quad(f2correção,0,np.inf) #integração
            #print(integral1)
            #print(integral1, T)
            B3corr2=((-N_A*10**-24)*(corr2i3))/(48*(irk**3)*(T**3)) #coef do virial da segunda correção
            
        
        # d1 = decimal.Decimal(B1)
        
            corr3i3, err3i3 =quad(f3correção,0,np.inf) #integração
            #print(integral1)
            #print(integral1, T)
            B3corr3=((-N_A*10**-24)*(corr3i3))/(48*(irk**3)*(T**3)) #coef do virial da terceira correção

        
            corr4i3, err4i3 =quad(f4correção,0,np.inf) #integração
            #print(integral1)
            #print(integral1, T)
            B3corr4=((-N_A*10**-24*(h**4))*(corr4i3))/(1920*(irk**4)*(T**4)*(mi**2)) #coef do virial da quarta correção
            
            
            Btotal3= (B3 + B3corr1 + B3corr2 + B3corr3 + B3corr4)  #soma de cada um dos coef encontrados resultando no coeficiente total
            
            #d1 = decimal.Decimal(B1)
            
            Tstr.append(T) #temperaturas usadas são colocadas na lista
            
            Bmonte.append(Bmonte3)
            Bstr.append(Btotal3) #Coef total achado é alocado na lista
            print('B=',Btotal3, 'para T=',T)
            

        r = np.linspace(0, 10, 100)


        integral_inf, err_inf =quad(lim_inf,0,np.inf) #integração
            #print(integral, err3)
        B_inf=(2*np.pi)*(integral_inf+(r0**3)/3) #coef do virial principal
        print('B_inf=', B_inf)

        integral_sup, err_sup =quad(lim_sup,0,np.inf) #integração
            #print(integral, err3)
        B_sup=(2*np.pi)*(integral_sup+(r0**3)/3) #coef do virial principal
        print('B_sup=', B_sup)
        #CRIANDO O GRÁFICO


        # Create two subplots and unpack the output array immediately
        f, (ax1, ax2) = plt.subplots(1, 2, constrained_layout=True)


        ax1.plot(Tstr, Bstr, label = 'dados calculados')
        ax1.plot(T_ref,B_ref,'.', label = 'dados de referência')
        ax1.plot(Tstr, Bmonte, label = 'dados monte carlo') #SEM CORREÇÕES QUANTICAS

        ax1.set_title('(a)')
        ax2.set_title('(b)')
        ax2.plot(r, U(r))

        ax1.set_xlim([-10, 500]) #1.2 B 
        ax1.set_ylim([-300, 35]) #inferior : B(T=10) * 1.1 superior: B(T=500)* 1.2
        ax2.set_xlim([0.3, 6])
        ax2.set_ylim([-1.5*De, 2*De]) #limite inferior De*1.5 #limite superior 3*De

        ax1.set_ylabel(r'$B(T)[cm^3/mol]$')
        ax1.set_xlabel(r'$Temperature[K]$', labelpad=1)
        ax2.set_ylabel(r'$U(r)[kcal/mol]$', labelpad=1)
        ax2.set_xlabel(r'$r [\AA]$',labelpad=1)
        '''
        ax_image = fig.add_subplot(Tstr, Bstr)
        ax_image.set_title('Imagem original')
        ax_image.imshow(image, cmap='gray')

        '''
        #plt.subplot(Tstr, Bstr)
        #plt.scatter(Tstr, Bstr)
        #plt.title(f'Gráficos (a) do , {name}. You are {age}.')

        ax1.legend()

        plt.show()