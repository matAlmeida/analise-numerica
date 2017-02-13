# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 21:49:40 2017

@author: evalero
"""

import numpy as np
import math as mt
import pylab
from edoNumericas import *

#Determinar y(t tal que)
#dy/dt = fde_(t,y) = -(y + 1)(y + 3)
def fde_(t,y):
    f = -(y + 1)*(y + 3)
    return f
#No intervalo 0 <= t <= 2.0
t_0 = 0.0
t_F = 2.0
#Com a cond. inicial y(t=0)=0,5
y_0 = 0.5
# A Solução analítica deste problema e
# y(t) = -3+2(1+e^-2t)^-1
def yde_(t):
    y_ = -3 + 2*(1 + mt.exp(-2 * t))**(-1)
    return y_

def dfdt_(t,y):
    df = y**3 + 6*(y**2) + 11*y + 6
    return df

def d2fdt2_(t,y):
    df = -3*(y**4) - 24*(y**3) - 68*(y**2) - 80*y - 33
    return df

def d3fdt3_(t,y):
    df = 3*(y**5) + 30*(y**4) + 115*(y**3) + 210*(y**2) + 182*y + 60
    return df


#Utilizando Taylor de 4ta grau com N = 10
N = 10
(t,w) = taylor4EDO(fde_, dfdt_, d2fdt2_, d3fdt3_, t_0, t_F, y_0, N)

pylab.figure(5)
pylab.subplot(211)
pylab.plot(t,w,'go--', label="Taylor 4ta")
pylab.title("Método de Runge-Kutta de 4ta ordem")
pylab.xlabel("t")
pylab.ylabel("y")

y = np.ones(t.size)
for i in range(0,t.size):
    y[i] = yde_(t[i])

print ("Método de Taylor 4ta Ordem com N=10")
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'go--', label="Taylor 4ta")

pylab.xlabel("t")
pylab.ylabel("Erro absoluto")

pylab.show()

#Utilizando RK de 3ra ordem com N = 10
N = 10
(t,w) = heunEDO3(fde_, t_0, t_F, y_0, N)

pylab.subplot(211)
pylab.plot(t,w,'yo--', label="RK 3ra")

print ("Método de Runge-Kutta 3ra Ordem com N=10")
errorAbs(t,y,w)
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'yo--', label="RK 3ra")

#Utilizando RK de 4ta ordem com N = 10
N = 10
(t,w) = rk4taOrdem(fde_, t_0, t_F, y_0, N)

pylab.subplot(211)
pylab.plot(t,w,'bo--', label="RK 4ta")

print ("Método de Runge-Kutta 4ta Ordem com N=10")
errorAbs(t,y,w)
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'bo--', label="RK 4ta")
pylab.legend(loc=2)

t = np.linspace(t_0, t_F, N*100)
y = np.ones(t.size)
for i in range(0,t.size):
    y[i] = yde_(t[i])

pylab.subplot(211)
pylab.plot(t,y,'r', label="y(t)")
pylab.legend(loc=2)
