# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 08:13:54 2017

@author: gabriel
"""

import numpy as np
import math as mt
import pylab
from edoNumericas import *

def fde_(t, y, j):
    g = 10.0
    L = 1.0
    return -(g/L) * mt.sin(y_0[j])

t_0 = 0.0
t_F = 10.0

global y_0
y_0 = np.zeros(2)
y_0[0] = 0
y_0[1] = mt.pi/4.0

N = 100
(t, w) = rk4taOrdemPlus(fde_, t_0, t_F, y_0, N)
#(t, w) = eulerEDOs(fde_, t_0, t_F, y_0, N)

pylab.figure(1)
pylab.subplot(211)
pylab.plot(t, w[0], 'go--', label="\teta")

pylab.subplot(212)
pylab.plot(t, w[1], 'go--', label="\alpha")
pylab.show()
