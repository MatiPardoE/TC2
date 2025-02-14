#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TS10 ej 1

@author: mariano
"""

import sympy as sp
import splane as tc2
from schemdraw import Drawing
from schemdraw.elements import  Resistor


# Resolución simbólica

s = sp.symbols('s ', complex=True)

# Sea la siguiente función de excitación
ZZ = (s**2+3)*(s**2+1)/(s*(s**2+2))

# a) ZZ según Foster derivación

# Implementaremos Imm mediante Foster
k0, koo, ki = tc2.foster(ZZ)

tc2.dibujar_foster_serie(k0, koo, ki, z_exc = ZZ)
# b) ZZ según Cauer1 (removiendo en oo) 

koo, imm_cauer_oo, rem = tc2.cauer_LC(ZZ, remover_en_inf=True)

tc2.print_latex( r'$' + sp.latex(ZZ) + r'=' + sp.latex(imm_cauer_oo) + r'$' )

# Tratamos a nuestra función inmitancia como una Z
tc2.dibujar_cauer_LC(koo, z_exc = imm_cauer_oo)

# # b) ZZ según Cauer2 (removiendo en 0) 
k0, imm_cauer_0, rem = tc2.cauer_LC(ZZ, remover_en_inf=False)

tc2.print_latex( r'$' + sp.latex(ZZ) + r'=' + sp.latex(imm_cauer_0) + r'$' )

# # Tratamos a nuestra función inmitancia como una Z
tc2.dibujar_cauer_LC(k0, z_exc = imm_cauer_0)
    


