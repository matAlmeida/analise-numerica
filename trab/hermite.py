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
        p = [(x - X[j]) / (X[i] - X[j]) for j in range(n) if j != i]
        return reduce(op.mul, p)

    def dL(i):
        d = [1.0 / (X[i] - X[j]) for j in range(n) if j!= i]
        return reduce(op.add, d)

    def A(i):
        return (1 - 2 * (x - X[i]) * dL(i)) * (L(i) ** 2)

    def B(i):
        return (x - X[i]) * (L(i) ** 2)

    assert(len(X) != 0 and len(X) == len(Y)), 'Quantidade de valores em X e Y diferentes'
    n = len(X)
    h1 = sum(A(i) * Y[i] for i in range(n))
    h2 = sum(B(i) * dY[i] for i in range(n))
    return h1 + h2

###########################################################################################

x = [1.3, 1.6, 1.9]
f_x = [0.6200860, 0.4554022, 0.2818186]
d_f_x= [-0.5220232, -0.5698959, -0.5811571]

print(H(1.5, x, f_x, d_f_x))
