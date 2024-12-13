import sympy as sp
import lib

#  Needed symbols
t, r, theta, phi, m = sp.symbols('t r theta phi m')
coords = [t, r, theta, phi]

m = sp.symbols('m')
V = sp.Function('V')(r)
U = sp.Function('U')(r)

# Metric tensor
metric = sp.Matrix([
    [U, 0, 0, 0],
    [0, -V, 0, 0],
    [0, 0, -r**2, 0],
    [0, 0, 0, -r**2 * sp.sin(theta)**2]
])

metric_inv, Gamma, Riemann, Ricci, Ricci_scalar, Einstein = lib.calculate_tensors(metric, coords)
# lib.pretty_print(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_scalar, Einstein)
lib.pretty_plt(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_scalar, Einstein)

# equations = []
# n = len(coords)
# for i in range(n):
#     for j in range(n):
#         if Einstein[i, j] != 0:
#             Einstein_component = Einstein[i, j]
#             equations.append(sp.Eq(Einstein_component, 0))
#
# solution = sp.solve(equations)
# print(solution)

# import matplotlib.pyplot as plt
#
# # Simple LaTeX rendering example in matplotlib
# def latex_demo():
#     latex_string = r'$\frac{d}{dx} e^x = e^x$'
#     plt.text(0.5, 0.5, latex_string, fontsize=15, ha='center', va='center')
#
#     plt.axis('off')
#     plt.tight_layout()
#     plt.show()
#
# # Run the demo
# latex_demo()
