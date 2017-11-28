############################################
#
# Implementación del método de Newton y
# primeras pruebas
#
# Autor: Francisco Luque Sánchez
#
############################################


def newton_method(f, df, x0, tolerance, iters):
    i = 0
    stop = False
    curr_x = x0
    while not stop:
        i += 1
        curr_val = f(curr_x)
        curr_x = curr_x - curr_val/df(x0)
        if abs(curr_val) < tolerance or i > iters:
            stop = True

    return curr_x


# Prueba del algoritmo con el polinomio
# f(x) = (x - 3)*(x - 2)*(x - 1) = x^3 - 6x^2 + 11x - 6
def f(x):
    return x**3 - 6*x**2 + 11*x - 6


def df(x):
    return 3*x**2 - 12*x + 11


for i in [0.5, 1.5, 2.5, 3.5]:
    x = newton_method(f, df, i, 10**-6, 10000000)
    print(x, f(x))
