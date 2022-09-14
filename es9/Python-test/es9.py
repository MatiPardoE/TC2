#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mariano

Descripción: Script para ejemplificar el uso de filtros digitales FIR e IIR 
estudiados en Teoría de Circuitos II. Se trabaja sobre una señal electrocardiográfica
registrada a 1 kHz, con diversos tipos de contaminación, que se buscan eliminar 
con los filtros diseñados. La plantilla de diseño se corresponde con un filtro
pasabanda con banda de paso de 3 a 25 Hz. Los detalles de la plantilla se pueden
ver en la implementación.

Algunos aspectos a considerar:
-----------------------------
    
    + El diseño de filtros FIR e IIR puede resolverse muy sencillamente con las rutinas incluidas
    en SciPy Signal. 
    
    + Sin embargo, algunas funciones son poco prácticas al momento de escribir estas rutinas.
    Las de filtros IIR funcionan correctamente, mientras que las de FIR dan bastante trabajo
    para cumplir con una especificación estricta, como la que mostramos en este ejemplo.
    
    + Los filtros diseñados, son efectivos para suprimir las interferencias, pero
    ineficaces para ser utilizados en la práctica, debido a la distorsión que
    introducen, especialmente a los tramos de señal no contaminados. Es decir, 
    son efectivos pero no inocuos.

"""

import scipy.signal as sig
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
from splane import plot_plantilla

def group_delay(ww, phase):
    
    groupDelay = -np.diff(phase)/np.diff(ww)
    
    return(np.append(groupDelay, groupDelay[-1]))


# Setup inline graphics
mpl.rcParams['figure.figsize'] = (10,10)

# para listar las variables que hay en el archivo
#io.whosmat('ecg.mat')
mat_struct = sio.loadmat('ecg.mat')

ecg_one_lead = mat_struct['ecg_lead']
ecg_one_lead = ecg_one_lead.flatten()
cant_muestras = len(ecg_one_lead)

fs = 1000 # Hz
nyq_frec = fs / 2


# filter design
ripple = 0.5 # dB
atenuacion = 40 # dB

ws1 = 1.0    #Hz
wp1 = 3.0    #Hz
wp2 = 25.0   #Hz
ws2 = 35.0   #Hz

frecs = np.array([0.0,         ws1,         wp1,     wp2,     ws2,         nyq_frec   ]) / nyq_frec
gains = np.array([-atenuacion, -atenuacion, -ripple, -ripple, -atenuacion, -atenuacion])
gains = 10**(gains/20)

bp_sos_cheby2 = sig.iirdesign(wp=np.array([wp1, wp2]) / nyq_frec, ws=np.array([ws1, ws2]) / nyq_frec, gpass=ripple, gstop=atenuacion, analog=False, ftype='cheby2', output='sos')

den = 1.0
demora = 1000
#*.csv diseñado en pyFDA
num_fir_lp = np.genfromtxt('LowPass500.csv' , delimiter = ',')
num_fir_hp = np.genfromtxt('HighPass1500.csv', delimiter = ',')
num_fir = sig.convolve(num_fir_lp, num_fir_hp)

w_rad  = np.append(np.logspace(-2, 0.8, 250), np.logspace(0.9, 1.6, 250) )
w_rad  = np.append(w_rad, np.linspace(40, nyq_frec, 500, endpoint=True) ) / nyq_frec * np.pi

_, h_cheby = sig.sosfreqz(bp_sos_cheby2, w_rad)
_, hh_fir = sig.freqz(num_fir, den, w_rad)


w = w_rad / np.pi * nyq_frec

plt.close('all')

plt.figure(1)
plt.clf()

plt.plot(w, 20*np.log10(np.abs(h_cheby)+1e-12), label='IIR-Cheby {:d}'.format(bp_sos_cheby2.shape[0]*2) )

plt.plot(w, 20 * np.log10(abs(hh_fir)), label='FIR-pyFDA {:d}'.format(num_fir.shape[0]))

plot_plantilla(filter_type = 'bandpass', fpass = frecs[[2, 3]]* nyq_frec, ripple = ripple , fstop = frecs[ [1, 4] ]* nyq_frec, attenuation = atenuacion, fs = fs)

plt.title('FIR diseñado por métodos directos')
plt.xlabel('Frequencia [Hz]')
plt.ylabel('Modulo [dB]')
plt.axis([0, 500, -60, 5 ]);

plt.grid()

axes_hdl = plt.gca()
axes_hdl.legend()
            
plt.figure(2)

phase_fir = np.angle(hh_fir)

plt.plot(w, phase_fir, label='FIR-pyFDA {:d}'.format(num_fir.shape[0]))    # Bode phase plot

plt.title('FIR diseñado retardo')
plt.xlabel('Frequencia [Hz]')
plt.ylabel('Fase [rad]')

axes_hdl = plt.gca()
axes_hdl.legend()

plt.show()

plt.figure(3)

gd_fir = group_delay(w_rad, phase_fir)

plt.plot(w[gd_fir>0], gd_fir[gd_fir>0], label='FIR-ls {:d}'.format(num_fir.shape[0]))    # Bode phase plot

plt.axis([0, 500, 0, 1.5*demora ])

plt.title('FIR diseñado retardo')
plt.xlabel('Frequencia [Hz]')
plt.ylabel('Retardo [s]')

axes_hdl = plt.gca()
axes_hdl.legend()

plt.show()

# #%% Ahora filtramos con cada filtro diseñado


# # FIltrado convencional

# # ECG_f_butt = sig.sosfilt(bp_sos_butter, ecg_one_lead)
# # ECG_f_cheb = sig.sosfilt(bp_sos_cheby, ecg_one_lead)
# # ECG_f_cauer = sig.sosfilt(bp_sos_cauer, ecg_one_lead)

# # ECG_f_ls = sig.lfilter(num_firls, den, ecg_one_lead)
# # ECG_f_remez = sig.lfilter(num_remez, den, ecg_one_lead)
# # ECG_f_win = sig.lfilter(num_win, den, ecg_one_lead)
# ECG_f_rl = sig.lfilter(num_rl, den_rl, ecg_one_lead)


# # # FIltrado bidireccional

# ECG_f_butt = sig.sosfiltfilt(bp_sos_butter, ecg_one_lead)
# ECG_f_cheb = sig.sosfiltfilt(bp_sos_cheby, ecg_one_lead)
# ECG_f_cauer = sig.sosfiltfilt(bp_sos_cauer, ecg_one_lead)

# ECG_f_ls = sig.filtfilt(num_firls, den, ecg_one_lead)
# ECG_f_remez = sig.filtfilt(num_remez, den, ecg_one_lead)
# ECG_f_win = sig.filtfilt(num_win, den, ecg_one_lead)
# # ECG_f_rl = sig.filtfilt(num_rl, den_rl, ecg_one_lead)





# regs_interes = ( 
#         np.array([5, 5.2]) *60*fs, # minutos a muestras
#         np.array([12, 12.4]) *60*fs, # minutos a muestras
#         np.array([15, 15.2]) *60*fs, # minutos a muestras
#         [4000, 5500], # muestras
#         [10e3, 11e3], # muestras
#         )

# for ii in regs_interes:
    
#     # intervalo limitado de 0 a cant_muestras
#     zoom_region = np.arange(np.max([0, ii[0]]), np.min([cant_muestras, ii[1]]), dtype='uint')
    
    
#     plt.figure()
#     plt.plot(zoom_region, ecg_one_lead[zoom_region], label='ECG', linewidth=2)
#     plt.plot(zoom_region, ECG_f_butt[zoom_region], label='Butter')
#     plt.plot(zoom_region, ECG_f_cheb[zoom_region], label='Cheby')
#     plt.plot(zoom_region, ECG_f_cauer[zoom_region], label='Cauer')
#     plt.plot(zoom_region, ECG_f_remez[zoom_region], label='Remez')
#     plt.plot(zoom_region, ECG_f_ls[zoom_region], label='LS')
#     plt.plot(zoom_region, ECG_f_win[zoom_region], label='Win')
    
#     # FIR con corrección de demora
#     # plt.plot(zoom_region, ECG_f_remez[zoom_region+demora], label='Remez')
#     # plt.plot(zoom_region, ECG_f_ls[zoom_region+demora], label='LS')
#     # plt.plot(zoom_region, ECG_f_win[zoom_region+demora], label='Win')
#     plt.plot(zoom_region, ECG_f_rl[zoom_region+demora_rl], label='Rick')
    
#     plt.title('ECG filtering example from ' + str(ii[0]) + ' to ' + str(ii[1]) )
#     plt.ylabel('Adimensional')
#     plt.xlabel('Muestras (#)')
    
#     axes_hdl = plt.gca()
#     axes_hdl.legend()
            
#     plt.show()



