'''
f) Filtro ecualizador de fase de 1ยบ orden. 

Verifique la transferencia que se obtendrรญa si a0=1, a1=-R, b0=R, b1=1 para 
R=-D/(D+2) y siendo D un valor de demora de -0,5 a 0,5 muestras (1/fs). 

En que valores de frecuencia este filtro obtendria un retardo de grupo acotado 
en un margen del 5% respecto a omega=0. Verificar que la demora obtenida es 
de 1+D muestras.
'''
# Inicialización e importación de módulos
# Módulos para Jupyter
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import scipy.signal as sig
import matplotlib as mpl
import matplotlib.pyplot as plt
from splane import GroupDelay

fig_sz_x = 10
fig_sz_y = 7
fig_dpi = 100 # dpi
fig_font_size = 16
mpl.rcParams['figure.figsize'] = (fig_sz_x,fig_sz_y)
plt.rcParams.update({'font.size':fig_font_size})

npoints = 500
offset = 0.1
fs = 2 
w_nyq = 2*np.pi*fs/2
d = [-0.5,0,0.5]

for current_d in d:

    r = (-current_d)/(current_d+2)
    print(r) #Para usar en pyFDA
    num_z = [r, 1]
    den_z = [1, r] #En la TF es (-a1) por eso r 
    eq_ph = sig.TransferFunction(num_z, den_z, dt=1/fs)
    #GroupDelay(eq_ph,filter_description="Ecualizador de fase D = {0}".format(current_d), digital = True, fs = fs)
    
    w, _, phase = eq_ph.bode(np.linspace(10**-2, w_nyq, npoints))
    phaseRad = phase * np.pi / 180.0
    groupDelay = -np.diff(phaseRad.reshape((npoints, 1)), axis = 0)/np.diff(w).reshape((npoints-1,1))
    groupDelay_correction = groupDelay * 2
    aux_hdl = plt.plot(w[1:] / w_nyq , groupDelay_correction, label="Ecualizador de fase D = {0}".format(current_d))# Bode phase plo
    
    margin_sup = groupDelay_correction[0][0] + groupDelay_correction[0][0] * (2.5/100)
    margin_inf = groupDelay_correction[0][0] - groupDelay_correction[0][0] * (2.5/100)
    
    for gd in range(len(groupDelay_correction)):
        if groupDelay_correction[gd][0] > margin_sup or groupDelay_correction[gd][0] < margin_inf:
            print('Omega Limite = F_nyq * {0} con D = {1}'.format(w[gd]/w_nyq,current_d))
            plt.annotate('Omega Lim 5% de demora = F_nyq * {0},'.format(w[gd]/w_nyq), xy=(w[gd]/w_nyq,groupDelay_correction[gd][0]), xytext=((w[gd]/w_nyq)+offset, groupDelay_correction[gd][0] + offset), arrowprops=dict(facecolor='black'))
            break
        if gd == len(groupDelay_correction)-1:
            print('Sin Omega Limite con D = {1}'.format(w[gd]/w_nyq, current_d))
            break
        
    
        
    
plt.grid(True)
plt.gca().set_xlim([0, 1])
plt.xlabel('Frecuencia normalizada a Nyq [#]')
plt.ylabel('Group Delay [sec]')
plt.title('Group delay')
plt.legend(loc='best')
axes_hdl = plt.gca()