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

from collections import namedtuple



st.title('Unit conversion calculator')

uploaded_file = st.file_uploader("upload a file")


input_unit = st.selectbox(
    'What is the input unit?',
    ('eV', 'meV', 'kcal-mol', 'cm^-1', 'hartree'))
output_unit = st.selectbox('What is the desired output unit?',
    ('eV', 'meV', 'kcal-mol', 'cm^-1', 'hartree'))

filename = st.text_input('write the name of the output file:')



def convert_mev():
    data = pd.read_csv(uploaded_file, sep="\s+", header=None)

    data.columns = ["R(Ang)", input_unit ]
    st.subheader('Input DataFrame')
    st.write(data)

    data[input_unit] = data[input_unit] - data[input_unit].iloc[-1]
    if output_unit == 'eV':
        data['eV'] = data[input_unit]/1000
        df_final = data[["R(Ang)", 'eV']]
        st.subheader('Output DataFrame')
        st.write(df_final)
        return df_final
    
    if output_unit == 'meV':
        df_final = data[["R(Ang)", 'meV']]
        st.subheader('Output DataFrame')
        st.write(df_final)
        return df_final
        
    elif output_unit == 'kcal-mol':
        data['eV'] = data[input_unit]/1000
        data['kcal-mol'] = data['eV']*23.0609
        df_final = data[["R(Ang)", 'kcal-mol']]
        st.subheader('Output DataFrame')
        st.write(df_final)
        return df_final

    elif output_unit == 'hartree':
        data['eV'] = data[input_unit]/1000
        data['hartree'] = data['eV'] * 0.0367502
        df_final = data[["R(Ang)", 'hartree']]
        st.subheader('Output DataFrame')
        st.write(df_final)
        return df_final
        
    elif output_unit == 'cm-1':
        data['eV'] = data[input_unit]/1000
        data['cm-1'] = data['eV'] * 8065.73
        df_final = data[["R(Ang)", 'cm-1']]
        st.subheader('Output DataFrame')
        st.write(df_final)
        return df_final


def convert_df(df):
    return df.to_csv(sep = " ", index = False,header = False).encode('utf-8')


 

if st.button('Convert'):
    start = time.time()
    if input_unit == 'meV':
        data = convert_mev()
    elif input_unit == 'eV':
        pass
    elif input_unit == 'cm-1':
        pass
    elif input_unit == 'kcal/mol':
        pass


    
    dat = convert_df(data)

    filename = filename + '.dat'

    st.download_button(
        "Press to Download your adjusted parameters",
        dat,
        filename,
        "text/csv",
        key='Download input file')



    st.write('convertion time {}'.format(time.strftime(
    "%H:%M:%S", time.gmtime(time.time()-start))))
    st.success('Calculations finished!', icon="✅")



