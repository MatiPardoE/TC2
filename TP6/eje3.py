# -*- coding: utf-8 -*-

import sympy as sp
import splane as tc2
from schemdraw import Drawing
from schemdraw.elements import  Resistor, Capacitor, Inductor
import numpy as np
from IPython import display

s = sp.symbols('s', complex=True)

ZZ = (2*s**2+1)/(s*(3*s**2+2))

# Restricción circuital: L1*C2 = 1/PI r/s
# remoción parcial en infinito de 1/YY

omega_L1C2 = round(np.sqrt(np.pi),5)

Y2, Yc1 = tc2.remover_polo_infinito(1/ZZ, omega_zero = omega_L1C2 )

print(omega_L1C2)
print(Yc1)

Y2_cool = Y2

for a in sp.preorder_traversal(Y2):
    if isinstance(a, sp.Float):
        Y2_cool = Y2_cool.subs(a, round(a, 5))

print(Y2_cool)

Z4, Zt2, L1, C2 = tc2.remover_polo_jw(1/Y2_cool, isImpedance = True, omega = omega_L1C2 )

display(Zt2)