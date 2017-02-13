# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 15:16:42 2017

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

#Utilizando RK de 2da ordem com N = 10
N = 10
(t,w) = rk2daOrdem(fde_, t_0, t_F, y_0, N)

pylab.figure(4)
pylab.subplot(211)
pylab.plot(t,w,'go--', label="RK 2da")
pylab.title("Métodos de Runge-Kutta de 2da ordem")
pylab.xlabel("t")
pylab.ylabel("y")

y = np.ones(t.size)
for i in range(0,t.size):
    y[i] = yde_(t[i])

print ("Método de Runge-Kutta 2da Ordem com N=10")    
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'go--', label="RK 2da")

pylab.xlabel("t")
pylab.ylabel("Erro absoluto")

#Utilizando RK (Euler Modificado) de 2da ordem com N = 10
N = 10 
(t,w) = eulerEDOMod(fde_, t_0, t_F, y_0, N) 

pylab.subplot(211)    
pylab.plot(t,w,'bo--', label="Euler Mod")
    
print ("Método de Euler Modificado com N=10")
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'bo--', label="Euler Mod")

#Utilizando RK (Heun) de 2da ordem com N = 10
N = 10 
(t,w) = heunEDO(fde_, t_0, t_F, y_0, N) 

pylab.subplot(211)    
pylab.plot(t,w,'yo--', label="Heun")
    
print ("Método de Heun com N=10")
e = errorAbs(t,y,w)

pylab.subplot(212)
pylab.plot(t,e,'yo--', label="Heun")
pylab.legend(loc=2)
    
t = np.linspace(t_0, t_F, N*100)
y = np.ones(t.size)
for i in range(0,t.size):
    y[i] = yde_(t[i])

pylab.subplot(211)     
pylab.plot(t,y,'r', label="y(t)")
pylab.legend(loc=2)
