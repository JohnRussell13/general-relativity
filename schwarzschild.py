import sympy as sp
import lib

#  Needed symbols
t, r, m, theta, phi = sp.symbols('t r m theta phi')
coords = [t, r, theta, phi]

# Metric tensor
metric = sp.Matrix([
    [-(1-2*m/r), 0, 0, 0],
    [0, (1-2*m/r)**(-1), 0, 0],
    [0, 0, r**2, 0],
    [0, 0, 0, r**2 * sp.sin(theta)**2]
])

metric_inv, Gamma, Riemann, Ricci, Ricci_curvature, Einstein = lib.calculate_tensors(metric, coords)
lib.pretty_print(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_curvature, Einstein)

# Manually went and proved that G_00 simplifies to 0

# equations = []
# n = len(coords)
# for i in range(n):
#     for j in range(n):
#         Einstein_component = Einstein[i, j]
#         
#         equations.append(sp.Eq(Einstein_component, 0))
#
# solution = sp.solve(equations)
#
# print(solution)
