# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

import numpy as np
import scipy.signal as sig
from splane import analyze_sys

fs = 1
num_z = [1,-1]
den_z = [1,0]

tf_digital = sig.TransferFunction(num_z,den_z,dt=1/fs)

analyze_sys([tf_digital],['digital'])