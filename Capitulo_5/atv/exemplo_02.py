# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 16:24:17 2017

@author: evalero
"""

import numpy as np
import math as mt
import pylab
from edoNumericas import *

#Exemplo 3 Capitulo 5 pag 233
#Determinar y(t tal que)
#dy/dt = fde_(t,y) = y - t^2 + 1
def fde_(t,y):
    f = y - t**2.0 + 1
    return f
#No intervalo 0 <= t <= 2.0
t_0 = 0.0
t_F = 2.0
#Com a cond. inicial y(t=0)=0,5
y_0 = 0.5
# A Solução analítica deste problema e
# y(t) = (t + 1)^2.0 - 0,5e^t
def yde_(t):
    y_ = (t + 1.0)**2.0 - 0.5*mt.exp(t)
    return y_

def dfdt_(t,y):
    df = y - t**2 - 2*t + 1
    return df

def dnfdtn_(t,y):
    df = y - t**2 - 2*t -1
    return df

#Utilizando Taylor de 2do grau com N = 10
N = 10
(t,w) = taylor2EDO(fde_, dfdt_, t_0, t_F, y_0, N)

pylab.figure(2)
pylab.subplot(211)
pylab.plot(t,w,'go--', label="Taylor 2da")
pylab.title("Método de Taylor")
pylab.xlabel("t")
pylab.ylabel("y")

y = np.ones(t.size)
for i in range(0,t.size):
    y[i] = yde_(t[i])

print ("Método de Taylor 2da Ordem com N=10")
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'go--', label="Taylor 2da")

pylab.xlabel("t")
pylab.ylabel("Erro absoluto")

pylab.show()

#Utilizando Taylor de 4to grau com N = 10
N = 10
(t,w) = taylor4EDO(fde_, dfdt_, dnfdtn_, dnfdtn_, t_0, t_F, y_0, N)

pylab.subplot(211)
pylab.plot(t,w,'bo--', label="Taylor 4ta")

print ("Método de Taylor 4ta Ordem com N=10")

errorAbs(t,y,w)
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'bo--', label="Taylor 4ta")
pylab.legend(loc=2)

t = np.linspace(t_0, t_F, N*100)
y = np.ones(t.size)
for i in range(0,t.size):
    y[i] = yde_(t[i])

pylab.subplot(211)
pylab.plot(t,y,'r', label="y(t)")
pylab.legend(loc=2)
