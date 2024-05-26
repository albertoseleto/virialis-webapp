from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import streamlit as st
from scipy.optimize import leastsq
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import sklearn.metrics  



st.title('Regression')
potential = st.selectbox(
    'What potential energy function do you want to use in the regression?',
    ('Rydberg Potential', 'Improved Leonard-Jonnes Potential'))

H_uploaded_file = st.file_uploader("upload a file containing the abinitio points for the H configuration")


guess = st.checkbox('I would like to make my own initial guess')


if potential == 'Rydberg Potential' :
    column_names = ['a1', 'a2', 'a3', 'a4', 'a5', 'De', 'Req', 'Eref', 'RMSE', 'correlation']
    row_names = ['H', 'X', 'Z', 'Ta', 'Tb', 'L']

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
        a1_guess =  3.5318
        a2_guess =  7.309
        a3_guess =  3.513
        a4_guess =  0.6
        a5_guess =  0.04
        De_guess =  0.0012
        Req_guess =  5.5
        Eref_guess =  -0.016

elif potential == 'Improved Leonard-Jonnes Potential' :
    column_names = ['alpha', 'beta', 'mp', 'De', 'Req','RMSE', 'correlation']
    row_names = ['H', 'X', 'Z', 'Ta', 'Tb', 'L']
    
    if guess:
        alpha_guess = st.text_input('write an initial guess for the value of alpha:(suggestion: 4)')
        beta_guess = st.text_input('write an initial guess for the value of beta:(suggestion: 8)')
        mp_guess = st.text_input('write an initial guess for the value of mp:(suggestion: 6)')
        De_guess = st.text_input('write an initial guess for the value of De:(suggestion: 52.86)')
        Req_guess = st.text_input('write an initial guess for the value of Req:(suggestion: 3.169)')

    else:
        alpha_guess = 0.5
        beta_guess = 0.53
        mp_guess = 3.3
        De_guess = 4000
        Req_guess = 0.8
try:
    a1_guess = float(a1_guess)
    a2_guess = float(a2_guess)
    a3_guess = float(a3_guess)
    a4_guess = float(a4_guess)
    a5_guess = float(a5_guess)
    De_guess = float(De_guess)
    Req_guess = float(Req_guess)
    Eref_guess = float(Eref_guess)
except:
    pass
try:
    alpha_guess = float(alpha_guess)
    beta_guess = float(beta_guess)
    mp_guess = float(mp_guess)  
    De_guess = float(De_guess)
    Req_guess = float(Req_guess)
except:
    pass

filename = st.text_input('write the name of the output file:')

## Função de otimização para regressão linear dos dados

if st.button('Calculate'):


    if potential == 'Rydberg Potential' :

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
            st.subheader('Fitting for H...')

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

            popt, pcov = curve_fit(rydberg, rH, UH, p0,maxfev=10000)
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
            



        def convert_df(df):
            return df.to_csv(sep = " ", index = False,header = False).encode('utf-8')

        dffinal = dfH




        st.subheader('Final table')
        df_info = dffinal

        dffinal = dffinal[['a1','a2','a3','a4','a5','De', 'Req','Eref']]
        st.write(dffinal)





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
        go.Scatter(x=rH, y=UH, name='H points', mode = 'markers')
                )
        fig1.add_trace(
        go.Scatter(x=rH, y=UfitH, name='H regression')
                    )
                

        fig1.add_trace(
        go.Scatter(x=rH, y=iniciais(rH), name='initial guess', mode = 'markers')
                    )



        fig1.update_layout(
            title_text="Graph of Simple Energies per Distance",
            yaxis_title='U[meV]',  # Replace 'X Axis Label' with your desired X-axis label
            xaxis_title=r"$r[A]$"   # Replace 'Y Axis Label' with your desired Y-axis label
                    )


    
        st.plotly_chart(fig1, use_container_width=True)

    elif potential == 'Improved Leonard-Jonnes Potential' :

        if alpha_guess is not None:
                    
            #st.write('a1',a1_guess, 'a2',a2_guess, 'a3',a3_guess, 'a4',a4_guess, 'a5', a5_guess, 'De',De_guess, 'Req',Req_guess, 'Eref',Eref_guess)

            def ILJ(r,alpha, beta,mp, De, Req):
                n = beta + alpha * (r / Req) ** 2
                return De * ((mp/(n - mp) * (Req/r) ** n) - (n/(n - mp) * (Req/r) ** mp))
                
        if H_uploaded_file is not None:

            data = pd.read_csv(H_uploaded_file, sep="\s+", header=None)
            data.columns = ["r", "U sem correc"] 
            st.subheader('Fitting for H')

            #De_guess = data["U sem correc"].min()
            #row = data[data["U sem correc"] == De_guess]
            #Req_guess = row['r'].iloc[0]

            #st.write("Minimum value from Column2:", De_guess)
            #st.write("Corresponding value from Column1:", Req_guess)

            values = [alpha_guess,beta_guess,mp_guess,De_guess,Req_guess]
            values = [float(v) for v in values]

            p0 = tuple(values)

            def iniciais(r):
                n = beta_guess + alpha_guess * (r / Req_guess) ** 2
                return De_guess * ((mp_guess/(n - mp_guess) * (Req_guess/r) ** n) - (n/(n - mp_guess) * (Req_guess/r) ** mp_guess))


            data['r'] = data['r']
            #st.write(data)

            rH = data['r']
            UH = data['U sem correc']

            popt, pcov = curve_fit(ILJ, rH, UH, p0,maxfev=5000)
            alpha,beta,mp, De, Req = popt
            yfitted = ILJ(rH,alpha,beta,mp, De, Req)



            
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

            UfitH=ILJ(rH,alpha, beta,mp, De, Req)

            st.write(dfH)

            #media = UH.mean()
            #dp = UH.std()
            #var = UH.var()

        def convert_df(df):
            return df.to_csv(sep = " ", index = False,header = False).encode('utf-8')

        dffinal = dfH


        st.subheader('Final table')

        df_info = dffinal

        dffinal = dffinal[['alpha', 'beta', 'mp', 'De', 'Req']]
        st.write(dffinal)

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
        go.Scatter(x=rH, y=UH, name='H points', mode = 'markers')
                )
        fig1.add_trace(
        go.Scatter(x=rH, y=UfitH, name='H regression')
                    )
        fig1.add_trace(
        go.Scatter(x=rH, y=iniciais(rH), name='initial guess', mode = 'markers')
                    )
        fig1.update_layout(
            title_text="Graph of Simple Energies per Distance",
            yaxis_title='U[eV]',  # Replace 'X Axis Label' with your desired X-axis label
            xaxis_title=r"$r[A]$"   # Replace 'Y Axis Label' with your desired Y-axis label
            )

        st.plotly_chart(fig1, use_container_width=True)
