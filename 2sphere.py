import sympy as sp
import lib

#  Needed symbols
R, theta, phi = sp.symbols('R theta phi')
coords = [theta, phi]

# Metric tensor
metric = sp.Matrix([
    [R**2, 0],
    [0, R**2 * sp.sin(theta)**2]
])

lib.calculate_print(metric, coords)
