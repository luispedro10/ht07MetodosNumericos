import numpy as np
from CalculoFunciones import *

def TerceraDerivadaAtras(f, a, h):
    return (1/h**3) * (f(a) - 3*f(a-h) + 3*f(a-2*h) - f(a-3*h))

def TerceraDerivadaCentrada(f, a, h):
    return (1/(2*h**3)) * (f(a+2*h) - 2*f(a+h) + 2*f(a-h) - f(a-2*h))

def f(x):
    return np.log(x**2)

h1 = 0.01
h2 = 0.001
x_val = 2

print("Tercera derivada usando fórmula Hacia atrás con h1:", TerceraDerivadaAtras(f, x_val, h1))
print("Tercera derivada usando fórmula Centrada con h1:", TerceraDerivadaCentrada(f, x_val, h1))
print("Tercera derivada usando fórmula Hacia atrás con h2:", TerceraDerivadaAtras(f, x_val, h2))
print("Tercera derivada usando fórmula Centrada con h2:", TerceraDerivadaCentrada(f, x_val, h2))


def calcular_error(aprox1, aprox2):
    return abs(aprox1 - aprox2)

error_atras = calcular_error(TerceraDerivadaAtras(f, x_val, h1), TerceraDerivadaAtras(f, x_val, h2))
error_centrada = calcular_error(TerceraDerivadaCentrada(f, x_val, h1), TerceraDerivadaCentrada(f, x_val, h2))

p_atras = np.log(error_atras / error_centrada) / np.log(h1 / h2)

print("Orden de convergencia para la fórmula Hacia atrás:", p_atras)
print("Orden de convergencia para la fórmula Centrada:", 3-p_atras)  # Dado que la suma de los órdenes es 3
