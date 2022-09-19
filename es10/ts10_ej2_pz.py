# -*- coding: utf-8 -*-

import sympy as sp
import splane as tc2
from schemdraw import Drawing
from schemdraw.elements import  Resistor, Capacitor, Inductor
import numpy as np

num = [3,0,7,0]
den = [1,0,7,0,10]

YY = tc2.TransferFunction(num,den)

tc2.pzmap(YY)