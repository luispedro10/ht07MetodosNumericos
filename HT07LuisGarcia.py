import numpy as np
from scipy.integrate import quad
from CalculoFunciones import *

# Definimos la función dada
def f(x):
    return np.cos(np.pi * x)**2

# Límites de integración
a = 0
b = 1/4

# Subintervalos
ms = np.array([4, 8, 16, 32])

# Aproximaciones con la Regla del Trapecio
print("Aproximaciones con la Regla del Trapecio Compuesta:")
for m in ms:
    Itrapecio = ReglaTrapecioCompuesta(f, a, b, m)
    print(f"m={m}, Itrapecio={Itrapecio:.6f}")

# Aproximaciones con la Regla de Simpson
print("\nAproximaciones con la Regla de Simpson Compuesta:")
for m in ms:
    Isimpson = ReglaSimpsonCompuesta(f, a, b, m)
    print(f"m={m}, Isimpson={Isimpson:.6f}")

# Valor exacto de la integral usando scipy
I_exacto, _ = quad(f, a, b)

# Mostrar el valor exacto
print(f"\nValor exacto de la integral: {I_exacto:.6f}")
