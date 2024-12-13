import sympy as sp

# Step 1: Define coordinates and the metric
t, r, theta, phi = sp.symbols('t r theta phi')  # Example: Spherical coordinates
coords = [t, r, theta, phi]

# Example metric: Schwarzschild metric in spherical coordinates
metric = sp.Matrix([
    [-(1 - 2 / r), 0, 0, 0],
    [0, 1 / (1 - 2 / r), 0, 0],
    [0, 0, r**2, 0],
    [0, 0, 0, r**2 * sp.sin(theta)**2]
])

# Step 2: Compute the inverse of the metric
metric_inv = metric.inv()

# Step 3: Calculate Christoffel symbols of the second kind
n = len(coords)
Gamma = sp.MutableDenseNDimArray.zeros(n, n, n)

for k in range(n):
    for i in range(n):
        for j in range(n):
            Gamma[k, i, j] = 1 / 2 * sum(
                metric_inv[k, l] * (sp.diff(metric[l, j], coords[i]) +
                                    sp.diff(metric[l, i], coords[j]) -
                                    sp.diff(metric[i, j], coords[l]))
                for l in range(n)
            )

# Step 4: Compute the Riemann tensor
Riemann = sp.MutableDenseNDimArray.zeros(n, n, n, n)

for l in range(n):
    for k in range(n):
        for i in range(n):
            for j in range(n):
                Riemann[l, k, i, j] = (sp.diff(Gamma[l, i, j], coords[k]) -
                                       sp.diff(Gamma[l, i, k], coords[j]) +
                                       sum(Gamma[m, i, j] * Gamma[l, k, m] -
                                           Gamma[m, i, k] * Gamma[l, j, m]
                                           for m in range(n)))

# Riemann tensor R^l_{kij}
for l in range(n):
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if Riemann[l, k, i, j] != 0:  # Show only non-zero components
                    print(f"Riemann[{l},{k},{i},{j}] = {Riemann[l, k, i, j]}")
