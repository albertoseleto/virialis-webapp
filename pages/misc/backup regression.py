from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import streamlit as st
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import sklearn.metrics  




potential = st.selectbox(
    'What potential energy function do you want to use in the regression?',
    ('Rydberg Potential', 'Improved Leonard-Jonnes Potential'))

H_uploaded_file = st.file_uploader("upload a file containing the abinitio points for the H configuration")
X_uploaded_file = st.file_uploader("upload a file containing the abinitio points for the X configuration")
Z_uploaded_file = st.file_uploader("upload a file containing the abinitio points for the Z configuration")
Ta_uploaded_file = st.file_uploader("upload a file containing the abinitio points for the Ta configuration")
Tb_uploaded_file = st.file_uploader("upload a file containing the abinitio points for the Tb configuration")
L_uploaded_file = st.file_uploader("upload a file containing the abinitio points for the L configuration")

guess = st.checkbox('I would like to make my own initial guess')
if guess:
    a1_guess = st.text_input('write an initial guess for the value of a1:(suggestion: 1.86)')
    a2_guess = st.text_input('write an initial guess for the value of a2:(suggestion: -1.5)')
    a3_guess = st.text_input('write an initial guess for the value of a3:(suggestion: 0.99)')
    a4_guess = st.text_input('write an initial guess for the value of a4:(suggestion: -0.277)')
    a5_guess = st.text_input('write an initial guess for the value of a5:(suggestion: 0.0371)')
    De_guess = st.text_input('write an initial guess for the value of De:(suggestion: 52.86)')
    Req_guess = st.text_input('write an initial guess for the value of Req:(suggestion: 3.169)')
    Eref_guess = st.text_input('write an initial guess for the value of Eref:(suggestion: 0.016)')

else:
    a1_guess =  1.86
    a2_guess =  -1.5
    a3_guess =  0.99
    a4_guess =  -0.277
    a5_guess =  0.0371
    De_guess =  52.86
    Req_guess =  3.169
    Eref_guess =  0.016

filename = st.text_input('write the name of the output file:')

## Função de otimização para regressão linear dos dados


column_names = ['a1', 'a2', 'a3', 'a4', 'a5', 'De', 'Req', 'Eref', 'RMSE', 'correlation']
row_names = ['H', 'X', 'Z', 'Ta', 'Tb', 'L']

if st.button('Calculate'):

    a1_guess = float(a1_guess)
    a2_guess = float(a2_guess)
    a3_guess = float(a3_guess)
    a4_guess = float(a4_guess)
    a5_guess = float(a5_guess)
    De_guess = float(De_guess)
    Req_guess = float(Req_guess)
    Eref_guess = float(Eref_guess)
    if a1_guess is not None:
                
        #st.write('a1',a1_guess, 'a2',a2_guess, 'a3',a3_guess, 'a4',a4_guess, 'a5', a5_guess, 'De',De_guess, 'Req',Req_guess, 'Eref',Eref_guess)
    

        
        def rydberg(r,a1,a2,a3,a4,a5,De, Req,Eref):
            y = (r-Req)
            U = -De*(1 + a1*y + a2*y**2 + a3*y**3 +
                        a4*y**4 + a5*y**5) * np.exp(-a1*y) + Eref
            return U
           
    if H_uploaded_file is not None:

        data = pd.read_csv(H_uploaded_file, sep="\s+", header=None)
        data.columns = ["r", "U sem correc"] 
        st.subheader('Fitting for H')

        #De_guess = data["U sem correc"].min()
        #row = data[data["U sem correc"] == De_guess]
        #Req_guess = row['r'].iloc[0]

        #st.write("Minimum value from Column2:", De_guess)
        #st.write("Corresponding value from Column1:", Req_guess)

        values = [a1_guess,a2_guess,a3_guess,a4_guess, a5_guess,De_guess,Req_guess,Eref_guess]
        values = [float(v) for v in values]

        p0 = tuple(values)

        def iniciais(r):
            y = (r-(Req_guess))/(Req_guess)
            U = -De_guess*(1 + a1_guess*y + a2_guess*y**2 + a3_guess*y**3 +
                        (a4_guess)*y**4 + a1_guess*y**5) * np.exp(-a1_guess*y) + Eref_guess
            return U


        data['r'] = data['r']
        #st.write(data)

        rH = data['r']
        UH = data['U sem correc']

        popt, pcov = curve_fit(rydberg, rH, UH, p0,maxfev=5000)
        a1,a2,a3,a4,a5,De, Req,Eref = popt
        yfitted = rydberg(rH,a1,a2,a3,a4,a5,De, Req,Eref)



        
        dfH = pd.DataFrame(popt)
        dfH = dfH.T

        mse = sklearn.metrics.mean_squared_error(UH, yfitted)  
        rmse = np.sqrt(mse)  
        dfH['RMSE'] = rmse
        corr = UH.corr(rH)
        dfH['correlation'] = corr
        dfH.columns = column_names

        row_H = ['H']
        dfH.index =row_H

        UfitH=rydberg(rH,a1,a2,a3,a4,a5,De, Req,Eref)

        st.write(dfH)

        #media = UH.mean()
        #dp = UH.std()
        #var = UH.var()
        

    if X_uploaded_file is not None:

        data = pd.read_csv(X_uploaded_file, sep="\s+", header=None)
        data.columns = ["r", "U sem correc"] #, "U com correc", "nothing", "nada","nitch"
        st.subheader('Fitting for X')


        #data = data.drop(["nothing", "nada","nitch"], axis = 1)
        data['r'] = data['r']
        #st.write(data)

        rX = data['r']
        UX = data['U sem correc']
        #p0 = (1.863598546970193,-1.529558156764738, 0.9953471906470752,-0.2776074439976977,0.03718004848996886,52.861310909897924,3.1691244646449928,0.0160366742839507)
        ## Otimização da curva


        popt, pcov = curve_fit(rydberg, rX, UX, p0,maxfev=5000)
        a1,a2,a3,a4,a5,De, Req,Eref = popt
        yfitted = rydberg(rX,a1,a2,a3,a4,a5,De, Req,Eref)

        dfX = pd.DataFrame(popt)
        dfX = dfX.T
        mse = sklearn.metrics.mean_squared_error(UX, yfitted)  
        rmse = np.sqrt(mse)  
        dfX['RMSE'] = rmse
        corr = UX.corr(rX)
        dfX['correlation'] = corr
        dfX.columns = column_names
        row_X = ['X']
        dfX.index =row_X

        UfitX=rydberg(rX,a1,a2,a3,a4,a5,De, Req,Eref)

        st.write(dfX)

    if Z_uploaded_file is not None:

        data = pd.read_csv(Z_uploaded_file, sep="\s+", header=None)
        data.columns = ["r", "U sem correc"] #, "U com correc", "nothing", "nada","nitch"
        st.subheader('Fitting for Z')


        #data = data.drop(["nothing", "nada","nitch"], axis = 1)
        data['r'] = data['r']
        #st.write(data)

        rZ = data['r']
        UZ = data['U sem correc']
        #p0 = (1.863598546970193,-1.529558156764738, 0.9953471906470752,-0.2776074439976977,0.03718004848996886,52.861310909897924,3.1691244646449928,0.0160366742839507)
        ## Otimização da curva


        popt, pcov = curve_fit(rydberg, rZ, UZ, p0,maxfev=5000)
        a1,a2,a3,a4,a5,De, Req,Eref = popt
        yfitted = rydberg(rZ,a1,a2,a3,a4,a5,De, Req,Eref)

        dfZ = pd.DataFrame(popt)
        dfZ = dfZ.T
        mse = sklearn.metrics.mean_squared_error(UZ, yfitted)  
        rmse = np.sqrt(mse)  
        dfZ['RMSE'] = rmse
        corr = UZ.corr(rZ)
        dfZ['correlation'] = corr
        dfZ.columns = column_names
        row_Z = ['Z']
        dfZ.index =row_Z

        UfitZ=rydberg(rZ,a1,a2,a3,a4,a5,De, Req,Eref)

        st.write(dfZ)

    if Ta_uploaded_file is not None:

        data = pd.read_csv(Ta_uploaded_file, sep="\s+", header=None)
        data.columns = ["r", "U sem correc"] #, "U com correc", "nothing", "nada","nitch"
        st.subheader('Fitting for Ta')


        #data = data.drop(["nothing", "nada","nitch"], axis = 1)
        data['r'] = data['r']
        #st.write(data)

        rTa = data['r']
        UTa = data['U sem correc']
        #p0 = (1.863598546970193,-1.529558156764738, 0.9953471906470752,-0.2776074439976977,0.03718004848996886,52.861310909897924,3.1691244646449928,0.0160366742839507)
        ## Otimização da curva


        popt, pcov = curve_fit(rydberg, rTa, UTa, p0,maxfev=5000)
        a1,a2,a3,a4,a5,De, Req,Eref = popt
        yfitted = rydberg(rTa,a1,a2,a3,a4,a5,De, Req,Eref)

        dfTa = pd.DataFrame(popt)
        dfTa = dfTa.T
        mse = sklearn.metrics.mean_squared_error(UTa, yfitted)  
        rmse = np.sqrt(mse)  
        dfTa['RMSE'] = rmse
        corr = UTa.corr(rTa)
        dfTa['correlation'] = corr
        dfTa.columns = column_names
        row_Ta = ['Ta']
        dfTa.index =row_Ta

        UfitTa=rydberg(rTa,a1,a2,a3,a4,a5,De, Req,Eref)

        st.write(dfTa)

    if Tb_uploaded_file is not None:

        data = pd.read_csv(Tb_uploaded_file, sep="\s+", header=None)
        data.columns = ["r", "U sem correc"] #, "U com correc", "nothing", "nada","nitch"
        st.subheader('Fitting for Tb')


        #data = data.drop(["nothing", "nada","nitch"], axis = 1)
        data['r'] = data['r']
        #st.write(data)

        rTb = data['r']
        UTb = data['U sem correc']
        #p0 = (1.863598546970193,-1.529558156764738, 0.9953471906470752,-0.2776074439976977,0.03718004848996886,52.861310909897924,3.1691244646449928,0.0160366742839507)
        ## Otimização da curva


        popt, pcov = curve_fit(rydberg, rTb, UTb, p0,maxfev=5000)
        a1,a2,a3,a4,a5,De, Req,Eref = popt
        yfitted = rydberg(rTb,a1,a2,a3,a4,a5,De, Req,Eref)

        dfTb = pd.DataFrame(popt)
        dfTb = dfTb.T
        mse = sklearn.metrics.mean_squared_error(UTb, yfitted)  
        rmse = np.sqrt(mse)  
        dfTb['RMSE'] = rmse
        corr = UTb.corr(rTb)
        dfTb['correlation'] = corr
        dfTb.columns = column_names
        row_Tb = ['Tb']
        dfTb.index =row_Tb

        UfitTb=rydberg(rTb,a1,a2,a3,a4,a5,De, Req,Eref)

        st.write(dfTb)


    if L_uploaded_file is not None:

        data = pd.read_csv(L_uploaded_file, sep="\s+", header=None)
        data.columns = ["r", "U sem correc"] #, "U com correc", "nothing", "nada","nitch"
        st.subheader('Fitting for L')


        #data = data.drop(["nothing", "nada","nitch"], axis = 1)
        data['r'] = data['r']
        #st.write(data)

        rL = data['r']
        UL = data['U sem correc']
        #p0 = (1.863598546970193,-1.529558156764738, 0.9953471906470752,-0.2776074439976977,0.03718004848996886,52.861310909897924,3.1691244646449928,0.0160366742839507)
        ## Otimização da curva


        popt, pcov = curve_fit(rydberg, rL, UL, p0,maxfev=5000)
        a1,a2,a3,a4,a5,De, Req,Eref = popt
        yfitted = rydberg(rL,a1,a2,a3,a4,a5,De, Req,Eref)




        dfL = pd.DataFrame(popt)
        dfL = dfL.T
        mse = sklearn.metrics.mean_squared_error(UL, yfitted)  
        rmse = np.sqrt(mse)  
        corr = UL.corr(rL)
        dfL['RMSE'] = rmse
        dfL['correlation'] = corr
        dfL.columns = column_names
        row_L = ['L']
        dfL.index =row_L

        UfitL=rydberg(rL,a1,a2,a3,a4,a5,De, Req,Eref)

        st.write(dfL)

    def convert_df(df):
        return df.to_csv(sep = " ", index = False,header = False).encode('utf-8')

    dffinal = dfH.merge(dfX, how = 'outer').merge(dfZ, how = 'outer').merge(dfTa, how = 'outer').merge(dfTb, how = 'outer').merge(dfL, how = 'outer')



    dffinal.columns = column_names
    dffinal.index = row_names


    st.subheader('Final table')
    st.write(dffinal)

    df_info = dffinal




    dat = convert_df(dffinal)

    filename = filename + '.dat'
    st.download_button(
    "Press to Download your adjusted parameters",
    dat,
    filename,
    "text/csv",
    key='Download input file')
    st.write('OBS: this is the input file you need to use virialis')
    st.write('OBS: downloading this file will reload the whole page!')






    df_info = df_info.reset_index()
    df_info = df_info.rename(columns = {'index': 'configuration'})
    df_info = pd.DataFrame(np.vstack([df_info.columns, df_info]))

    st.subheader('Information table')

    st.write(df_info)

    filename2 = filename + '_info.dat'
    st.download_button(
    "Press to Download your information file",
    dat,
    filename2,
    "text/csv",
    key='Download')
    st.write('OBS: this file contains the input information of virialis, as well as information to help you remember what each value means')
    st.write('OBS: downloading this file will reload the whole page!')

    fig1 = go.Figure()

    fig1.add_trace(
    go.Scatter(x=rH, y=UH, name='pontos H', mode = 'markers')
            )
    fig1.add_trace(
    go.Scatter(x=rH, y=UfitH, name='ajuste H')
                )
            


    fig1.add_trace(
    go.Scatter(x=rX, y=UX, name='pontos X', mode = 'markers')
            )
    fig1.add_trace(
    go.Scatter(x=rX, y=UfitX, name='ajuste X')
                )

    fig1.add_trace(
    go.Scatter(x=rZ, y=UZ, name='pontos Z', mode = 'markers')
            )
    fig1.add_trace(
    go.Scatter(x=rZ, y=UfitZ, name='ajuste Z')
                )
    fig1.add_trace(
    go.Scatter(x=rTa, y=UTa, name='pontos Ta', mode = 'markers')
            )
    fig1.add_trace(
    go.Scatter(x=rTa, y=UfitTa, name='ajuste Ta')
                )

    fig1.add_trace(
    go.Scatter(x=rTb, y=UTb, name='pontos Tb', mode = 'markers')
            )
    fig1.add_trace(
    go.Scatter(x=rTb, y=UfitTb, name='ajuste Tb')
                )

    fig1.add_trace(
    go.Scatter(x=rL, y=UL, name='pontos L', mode = 'markers')
            )
    fig1.add_trace(
    go.Scatter(x=rL, y=UfitL, name='ajuste L')
                )
    fig1.add_trace(
    go.Scatter(x=rH, y=iniciais(rH), name='initial guess', mode = 'markers')
                )



    st.plotly_chart(fig1, use_container_width=True)

