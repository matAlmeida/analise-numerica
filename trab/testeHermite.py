#from scipy import interpolate
#from numpy import polyder
from functools import reduce
import operator as op

def H(x, X, Y, dY):
    """
    x - Valor a ser calculado
    X - Vetor com valores de x
    Y - Vetor com valores de f(x)
    dY - Vetor com valores de f'(x)
    """

    def L(i):
        #return p[i] * (x ** i)
        p = [(x - X[i]) / (X[j] - X[i]) for j in range(n) if j != i]
        return reduce(op.mul, p)

    def dL(i):
        #return d[i-1] * (x ** (i-1))
        if i < n-1:
            return (Y[i+1] - Y[i]) / (X[i+1] - X[i])
        else:
            return (Y[i] - Y[i-1]) / (X[i] - X[i-1])

    def A(i):
        return (1 - 2 * (x - X[i]) * dL(i)) * (L(i) ** 2)

    def B(i):
        return (x - X[i]) * (L(i) ** 2)

    assert(len(X) != 0 and len(X) == len(Y)), 'Quantidade de valores em X e Y diferentes'
    n = len(X)
    #p = interpolate.lagrange(X, Y)
    #d = polyder(p)
    h1 = sum(A(i) * Y[i] for i in range(n))
    h2 = sum(B(i) * dY[i] for i in range(n))
    return h1 + h2

###########################################################################################

x = [1.3, 1.6, 1.9]
f_x = [0.6200860, 0.4554022, 0.2818186]
d_f_x= [-0.5220232, -0.5698959, -0.5811571]

#p = interpolate.lagrange(x, f_x)
#print(p)
#print(p(1.5))
#print(p[-1])
print(H(1.5, x, f_x, d_f_x))
