import sympy as sp
import lib

#  Needed symbols
R, t, theta, phi = sp.symbols('r t theta phi')
coords = [t, theta, phi]

# Metric tensor
metric = sp.Matrix([
    [-1, 0, 0],
    [0, R**2, 0],
    [0, 0, R**2 * sp.sin(theta)**2]
])

lib.calculate_print(metric, coords)
