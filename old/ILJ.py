import streamlit as st
import numpy as np
import altair as alt
import pandas as pd

import vegas
import random
import time
import math
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from statistics import mean
from numpy import arange
from pandas import read_csv
from matplotlib import pyplot



alfa = 5.4392	
beta = 7.63703	
mp = 5.67314	
Req = 3.21477	
De = 12.5264 

def U_ILJ(r):

    n = beta + alfa * (r/ Req) ** 2
    return De * ((mp/(n - mp) * (r/Req) ** n) - (n/(n - mp) * (r/Req) ** mp))

    
r = np.linspace(0, 10, 100)


plt.plot(r, U_ILJ(r), label = 'UH')
plt.ylim([-1000,1000])

plt.show()



