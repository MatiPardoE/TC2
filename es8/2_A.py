# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

import numpy as np
import scipy.signal as sig
from splane import analyze_sys

fc = 1e3
wo = 2*np.pi*fc
norm_w = wo
norm_f = fc
fsampling_1 = 100e3
fsampling_2 = 10e3

wo_norm = wo/norm_w
fsampling_1_norm = fsampling_1/norm_f
fsampling_2_norm = fsampling_2/norm_f

q_butter = 1/np.sqrt(2)

num = [wo_norm]
den = [1,wo_norm/q_butter,wo_norm**2]

tf_analogica = sig.TransferFunction(num,den)

numz_1, denz_1 = sig.bilinear(num, den, fs=fsampling_1_norm)
tf_digital_1 = sig.TransferFunction(numz_1,denz_1,dt=1/fsampling_1_norm)

numz_2, denz_2 = sig.bilinear(num, den, fs=fsampling_2_norm)
tf_digital_2 = sig.TransferFunction(numz_2,denz_2,dt=1/fsampling_2_norm)



analyze_sys([tf_analogica,tf_digital_1,tf_digital_2],['Analogica','Digital Fs = 100KHz','Digital Fs = 10KHz'])




