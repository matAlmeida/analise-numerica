# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 21:14:13 2017

@author: evalero
"""

import numpy as np
import math as mt
import pylab
from edoNumericas import *

#Exemplo 1 Capitulo 5 pag 225
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

#Utilizando N = 10
N = 10
(t,w) = eulerEDO(fde_, t_0, t_F, y_0, N) 

pylab.figure(1)
pylab.subplot(211)
pylab.plot(t,w,'go--', label="N=10")
pylab.title("Método de Euler")
pylab.xlabel("t")
pylab.ylabel("y")

y = np.ones(t.size)
for i in range(0,t.size):
    y[i] = yde_(t[i])

print ("Método de Euler com N=10")    
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'go--', label="N=10")

pylab.xlabel("t")
pylab.ylabel("Erro absoluto")

pylab.show()

#Utilizando agora N = 20
N = 20 
(t2,w2) = eulerEDO(fde_, t_0, t_F, y_0, N) 

pylab.subplot(211)    
pylab.plot(t2,w2,'bo--', label="N=20")

for i in range(0,w.size):
    w[i] = w2[2*i]

print ("Método de Euler com N=20")
errorAbs(t,y,w)
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'bo--', label="N=20")
pylab.legend(loc=2)
    
t = np.linspace(t_0, t_F, N*100)
y = np.ones(t.size)
for i in range(0,t.size):
    y[i] = yde_(t[i])

pylab.subplot(211)    
pylab.plot(t,y,'r', label="y(t)")
pylab.legend(loc=2)
