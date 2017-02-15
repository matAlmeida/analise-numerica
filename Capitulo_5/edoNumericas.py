# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 22:20:58 2017

@author: evalero
"""

import numpy as np
import math as mt

# Método de Euler
def eulerEDO(f, t0, tF, y0, N):
    h = (tF - t0)/N
    t = np.linspace(t0, tF, N+1)
    w = np.zeros(t.size);

    w[0] = y0
    for i in range(0,t.size-1):
        w[i+1] = w[i]+h*f(t[i],w[i])

    return (t,w)

def eulerEDOs(f, t0, tF, y0, N):
    h = (tF - t0)/N
    t = np.linspace(t0, tF, N+1)

    eSize = y0.size

    s = (t.size, y0.size)

    w = np.zeros(t.size);

    for j in range(eSize):
        w[0][j] = y0[j]

    for i in range(0,t.size-1):
        for j in range(eSize):
            w[i+1][j] = w[i][j]+h*f(t[i],w[i][j])

    return (t,w)

def taylor2EDO(f, dfdt, t0, tF, y0, N) :
    h = (tF - t0)/N
    t = np.linspace(t0, tF, N+1)
    w = np.zeros(t.size);

    w[0] = y0
    for i in range(0,t.size-1):
        w[i+1] = w[i] + h* (f(t[i],w[i]) + 0.5*h*dfdt(t[i],w[i]))

    return (t,w)

def taylor4EDO(f, dfdt, d2fdt2, d3fdt3, t0, tF, y0, N) :
    h = (tF - t0)/N
    t = np.linspace(t0, tF, N+1)
    w = np.zeros(t.size);

    w[0] = y0
    for i in range(0,t.size-1):
        w[i+1] = f(t[i],w[i]) + 0.5*h*dfdt(t[i],w[i])
        w[i+1] += ((h**2.0)*d2fdt2(t[i],w[i]))/6.0
        w[i+1] += ((h**3.0)*d3fdt3(t[i],w[i]))/24.0
        w[i+1] *= h
        w[i+1] += w[i]
    return (t,w)

def rk2daOrdem(f, t0, tF, y0, N):
    h = (tF - t0)/N
    hh = 0.5*h
    t = np.linspace(t0, tF, N+1)
    w = np.zeros(t.size);

    w[0] = y0
    for i in range(0,t.size-1):
        w[i+1] = w[i] + h*f(t[i] + hh, w[i] + hh*f(t[i],w[i]))

    return (t,w)

def eulerEDOMod(f, t0, tF, y0, N):
    h = (tF - t0)/N
    hh = 0.5*h
    t = np.linspace(t0, tF, N+1)
    w = np.zeros(t.size);

    w[0] = y0
    for i in range(0,t.size-1):
        w[i+1] = w[i] + hh*(f(t[i], w[i])+f(t[i+1], w[i] + h*f(t[i], w[i])))

    return (t,w)

def heunEDO(f, t0, tF, y0, N):
    h = (tF - t0)/N
    qh = 0.25*h
    dth = 2.0*h/3.0
    t = np.linspace(t0, tF, N+1)
    w = np.zeros(t.size);

    w[0] = y0
    for i in range(0,t.size-1):
        w[i+1] = w[i] + qh*(f(t[i], w[i]) + 3*f(t[i]+dth, w[i] + dth*f(t[i], w[i])))

    return (t,w)

def heunEDO3(f, t0, tF, y0, N):
    h = (tF - t0)/N
    qh = 0.25*h
    dth = 2.0*h/3.0
    uth = h/3.0
    t = np.linspace(t0, tF, N+1)
    w = np.zeros(t.size);

    w[0] = y0
    for i in range(0,t.size-1):
        w[i+1] = w[i] + qh*(f(t[i], w[i]) + 3*f(t[i]+dth, w[i] + dth*f(t[i]+uth, w[i]+uth*f(t[i],w[i]))))

    return (t,w)

def rk4taOrdem(f, t0, tF, y0, N):
    h = (tF - t0)/N
    hh = 0.5*h
    t = np.linspace(t0, tF, N+1)
    w = np.zeros(t.size);

    w[0] = y0
    for i in range(0,t.size-1):
        k1 = h*f(t[i], w[i])
        k2 = h*f(t[i]+hh, w[i]+0.5*k1)
        k3 = h*f(t[i]+hh, w[i]+0.5*k2)
        k4 = h*f(t[i+1], w[i]+k3)
        w[i+1] = w[i] + (k1 + 2*k2 + 2*k3 + k4)/6.0

    return (t,w)

def rk4taOrdemPlus(f, t0, tF, y0, N):
    h = (tF - t0)/N
    hh = 0.5*h
    t = np.linspace(t0, tF, N+1)
    tam = y0.size
    w = []

    for j in range(tam):
        w.append([])
        w[j].append(float(y0[j]))
        for i in range(0,t.size-1):
            k1 = h*f(t[i], w[j][i], j)
            k2 = h*f(t[i]+hh, w[j][i]+0.5*k1, j)
            k3 = h*f(t[i]+hh, w[j][i]+0.5*k2, j)
            k4 = h*f(t[i+1], w[j][i]+k3, j)
            k = w[j][i] + (k1 + 2*k2 + 2*k3 + k4)/6.0
            w[j].append(k)

    return (t,w)

def errorAbs(t,y,w):
    e = np.zeros(t.size);
    for i in range(0,y.size):
        e[i] = mt.fabs(y[i]-w[i])
        print("%f - %f, - %f - %f \n" % (t[i], y[i], w[i], e[i]))
    return e
