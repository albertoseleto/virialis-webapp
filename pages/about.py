import streamlit as st
import numpy as np
import time
st.header('About Viriális')

st.write('Viriális is a code that calculates the second virial coefficient.')
st.write('It is created by undergraduate University of São Paulo student **Alberto Seleto de Souza**, through a CNPq scholarship program in INPE. The coordenator is **Dr. Patricia Regina Pereira Barreto**.')


st.subheader('How the Code Works')


st.write(r"""To begin with, it is nice to review the idea behind calculating the Second Virial Coefficient, and where it comes from. 
When working with real gases, equations like the Clapeyron equation does not work to determine the gases' thermodynamics properties. That is why, it is necessary to use the Virial Equation, such as follows:
$$
\frac{P V_{\mathrm{m}}}{R T}=1+\frac{B(T)}{V_{\mathrm{m}}}+\frac{C(T)}{V_{\mathrm{m}}^{2}}+\frac{D(T)}{V_{\mathrm{m}}^{3}}+\cdots
$$
Where $P$ is pressure, $T$ is temperature, $R$ is the gas constant, $V_{\mathrm{m}}$ is the molar volume, and $B(T), C(T), D(T)$ are the second, third, and fourth virial coefficients.

For now and probably forever, I will focus only on the Second Virial Coefficient. It is given by:

$$
B(T)= \frac{N_{A}}{4} \int_{0}^{2 \pi} \int_{0}^{\pi} \sin \theta_{1} \int_{0}^{\pi} \sin \theta_{2}
 \int_{0}^{\infty}\left(1-\exp \left(-\frac{V}{k_{B} T}\right)\right) r^{2} d r d \theta_{1} d \theta_{2} d \phi
$$

Where $V$ is the potential energy which depends on the coordinates ($r$, $\theta_{1}$, $\theta_{2}$, $\phi$), $k_{b}$ is the Boltzmann constant, and $N_{A}$ avogadro's number.

Explain the harmonics expansion part????

Now, we need to decide which potential energy function to use. In Viriális, the first option is the Rydberg Potential:

$$
V(r, \gamma)=-D_{min}(\gamma)\left[1+\left(\sum_{i=1}^{5} a_{i}(\gamma)\left(r-r_{min}(\gamma)\right)^{i}\right)\right] e^{-a_{1}(\gamma)\left(r-r_{min}(\gamma)\right)}+E_{ref}(\gamma)
$$

Where $\gamma \equiv\left(\theta_{1}, \theta_{2}, \phi\right)$, $D_{min}(\gamma)$ is the minimum energy, $r_{min}$ is the distance at the minimum energy point, $a_{i}(\gamma)$ are the polinomial expansion coefficients., e $E_{ref}(\gamma)$ is the reference energy.

         
The other option available is the Improved Leonnard-Jonnes Potential, defined as:
         
$$
V(r, \gamma)=D_{\min }(\gamma)\left[\frac{m_{P}(\gamma)}{n(r, \gamma)-m_{P}(\gamma)}\left(\frac{r_{\min }(\gamma)}{r}\right)^{n(r, \gamma)}-\frac{n(r, \gamma)}{n(r, \gamma)-m_{P}(\gamma)}\left(\frac{r_{\min }(\gamma)}{r}\right)^{m_{P}(\gamma)}\right]
$$
         
Where 
$$
n(r, \gamma)=b(\gamma)+a(\gamma)\left(\frac{r}{r_{\min }(\gamma)}\right)^{2}
$$

         
""")


st.subheader('work in progress...')
