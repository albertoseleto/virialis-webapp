from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np

r = np.linspace(0, 10, 100)
def UH(r):
    return r**2
def UHH(r):
    return r**3

fig = make_subplots(rows=1, cols=2)

fig.add_trace(
    go.Scatter(x=r, y=UH(r)),
    row=1, col=1
)
fig.add_trace(
    go.Scatter(x=r, y=UHH(r)),
    row=1, col=1
)

fig.add_trace(
    go.Scatter(x=r, y=UH(r)),
    row=1, col=2
)

fig.update_layout(height=600, width=800, title_text="Side By Side Subplots")
fig.show()






               f, (ax1, ax2,ax3) = plt.subplots(1, 3, constrained_layout=True)
               ax1.plot(r, LC.UH(r), label = 'UH')
               ax1.plot(r, LC.UX(r), label = 'UX')
               ax1.plot(r, LC.UZ(r), label = 'UZ')
               ax1.plot(r, LC.UTa(r), label = 'UTa')
               ax1.plot(r, LC.UTb(r), label = 'UTb')
               ax1.plot(r, LC.UL(r), label = 'UL')
               ax1.set_title('(a)energias simples')
               #ax1.set_xlim([2, 8]) #1.2 B 
               #ax1.set_ylim([-40, 100]) #inferior : B(T=10) * 1.1 superior: B(T=500)* 1.2
               ax1.set_ylabel(r'$U(r)[eV]$')
               ax1.set_xlabel(r'$r[Å]$', labelpad=1)
               #ax2.set_xlim([2, 8]) #1.2 B 
               #ax2.set_ylim([-40, 100])
               #ax3.set_xlim([2, 10]) #1.2 B 
               #ax3.set_ylim([-1000, 1006])
               #ax4.set_xlim([2, 10]) #1.2 B 
               #ax4.set_ylim([-1000, 1006])



               ax2.plot(r, complexa.UM_000(r),  color='black',label = 'U000')
               ax2.plot(r, complexa.UM_202(r), color='green',label = 'U202')
               ax2.plot(r, complexa.UM_220(r), color='blue',label = 'U220')
               ax2.plot(r, complexa.UM_022(r), color='red',label = 'U022')
               ax2.plot(r, complexa.UM_222(r), color='cyan',label = 'U222')
               ax2.plot(r, complexa.UM_224(r), color='magenta',label = 'U224')
               ax2.set_title('(b)energias complexas')
               #ax2.set_xlim([2, 9]) #1.2 B 
               #ax2.set_ylim([-0.4, 2]) #inferior : B(T=10) * 1.1 superior: B(T=500)* 1.2
               ax2.set_ylabel(r'$U(r)[eV]$')
               ax2.set_xlabel(r'$r[Å]$', labelpad=1)

               ax3.plot(Tstr, B_clas, color='r',label = 'Bclássico')
               #ax3.plot(T_ref, B_ref,   color='b',label = 'ref') #Breferencia
               ax3.plot(T_state, B_virial_state_ref,   color='g',label = 'Bvirialstate')



               #ax3.plot(Tstr, B_main,   color='yellow',label = 'Bmain')
               ax3.plot(Tstr, B_plus_all_except_c2, color='blue',label = 'Btotal')

               ax3.set_title('(c) Segundo Coeficiente Virial(B) em função da Temperatura(T)')
               ax3.set_ylabel(r'$B[cm^3/mol]$')
               ax3.set_xlabel(r'$T[K]$', labelpad=1)
               '''
               ax4.plot(Tstr, B_c1, color='g',label = 'Bc1')
               ax4.plot(Tstr, B_c2,  color='red',label = 'Bc2')
               ax4.plot(Tstr, B_c3,   color='blue',label = 'Bc3')
               ax4.plot(Tstr, B_c4,   color='magenta',label = 'Bc4')

               #ax3.set_ylim([-100000, 100000]) #1.2 B 



               ax4.plot(Tstr, B_correcoes,   color='black',label = 'Bcorreções')



               ax_image = fig.add_subplot(Tstr, Bstr)
               ax_image.set_title('Imagem original')
               ax_image.imshow(image, cmap='gray')
               '''

               #plt.subplot(Tstr, Bstr)
               #plt.scatter(Tstr, Bstr)
               #plt.title(f'Gráficos (a) do , name. You are age.')

               ax1.legend()
               ax2.legend()
               ax3.legend()
               #ax4.legend()