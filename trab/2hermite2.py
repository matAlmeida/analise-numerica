# -*- coding: utf-8 -*-
#
# interpolation_hermite.py - Hermite Interpolation to Fit Curves.
#
#
# Autor: Pedro Garcia Freitas [sawp@sawp.com.br]
# License: Creative Commons
#      <http://creativecommons.org/licenses/by-nc-nd/2.5/br/>
#
# Dec 2010
#


def interpolation_hermite(x, x_j, f_j, diffF_j, n=0):
    """
    Return x evaluated in Hermite interpolation function.

    interpolation_hermite(x, x_j, f_j, diffF_j, n=0)

    INPUT:
      * x: evalue at this point
      * x_j: LIST, tabular points
      * f_j: LIST, tabular points (must be same size of x_j)
      * diffF_j: LIST, derivative points of x_j (f'(x_j))
      * n: polynomial degree (optional)

    return:
      * y(x): interpolated and evaluated x value

    Author: Pedro Garcia <sawp@sawp.com.br>
    see: http://www.sawp.com.br

    License: Creative Commons
             http://creativecommons.org/licenses/by-nc-nd/2.5/br/

    Dec 2010
    """

    if type(x_j) != type(f_j) or type(x_j) != type([]):
        print "Error: wrong type parameters"
        return float("NaN")

    if len(x_j) != len(f_j) or len(x_j) != len(diffF_j):
        print "Error: the tabular points must have same size"
        return float("NaN")

    if n <= 0:
        n = len(x_j)
    else:
        n = n + 1

    p = 0.0
    for j in xrange(n):
        # lagrange polinomial
        xj = x_j[j]
        ljn_num = 1.0
        ljn_den = 1.0
        for i in xrange(n):
            xi = x_j[i]
            if i != j:
                ljn_num *= (x - xi)
                ljn_den *= (xj - xi)
        ljn = ljn_num / ljn_den

        # hermite interpolation
        diff_ljn = 0
        for i in xrange(n):
            xi = x_j[i]
            if i != j:
                diff_ljn += 1.0 / (xj - xi)

        hjn = (1.0 - 2 * (x - xj) * diff_ljn) * (ljn * ljn)
        hjn_ = (x - xj) * (ljn * ljn)
        p += hjn * f_j[j] + hjn_ * diffF_j[j]
    return p

if __name__ == "__main__":
    J0 = [1.0]
    J0.append(0.7951976866)
    J0.append(0.2238907791)
    J0.append(-0.2600519549)
    J0.append(-0.3971498099)
    J0.append(-0.1775967713)
    J0.append(0.1506452573)
    J0.append(0.3000792705)
    J0.append(0.1716508071)
    J0.append(-0.0903336112)
    J0.append(-0.2459357645)

    # in this case, f == J1 (bessel 1th order)
    f = [0.0]
    f.append(0.4400505857)
    f.append(0.5767248078)
    f.append(0.3390589585)
    f.append(-0.0660433280)
    f.append(-0.3275791376)
    f.append(-0.2766838581)
    f.append(-0.0046828235)
    f.append(0.2346363469)
    f.append(0.2453117866)
    f.append(0.0434727462)

    J2 = [0.0]
    J2.append(0.1149034849)
    J2.append(0.3528340286)
    J2.append(0.4860912606)
    J2.append(0.3641281459)
    J2.append(0.0465651163)
    J2.append(-0.2428732100)
    J2.append(-0.3014172201)
    J2.append(-0.1129917204)
    J2.append(0.1448473415)
    J2.append(0.2546303137)

    termos = []
    for i in xrange(len(f)):
        termos.append(float(i))

    # in this problem, evaluate J1', where J1' = (J0 - J2)/2
    # so, diffF = J1'
    diffF = []
    for i in xrange(len(f)):
        k = (J0[i] - J2[i]) / 2.0
        diffF.append(k)

    # J
    #for i in xrange(1, 22):
    #    interp = interpolation_hermite(i * 0.5, termos, f, diffF)
    #    print "J(", i * 0.5, ") = ", interp
    x = [1.3, 1.6, 1.9]
    f_x = [0.6200860, 0.4554022, 0.2818186]
    d_f_x= [-0.5220232, -0.5698959, -0.5811571]
    p = interpolation_hermite(1.5, x, f_x, d_f_x)
    print p
