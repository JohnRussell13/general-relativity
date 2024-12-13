import sympy as sp

# Step 0: Needed symbols
r, theta, phi = sp.symbols('r theta phi')
coords = [r, theta, phi]

# Step 1: Metric tensor
metric = sp.Matrix([
    [1, 0, 0],
    [0, r**2, 0],
    [0, 0, r**2 * sp.sin(theta)**2]
])

# Step 2: Inverese Metric tensor
metric_inv = metric.inv()

# Step 3: Christoffel symbols
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

# Step 4: Riemann tensor
Riemann = sp.MutableDenseNDimArray.zeros(n, n, n, n)

for l in range(n):
    for k in range(n):
        for i in range(n):
            for j in range(n):
                Riemann[l, k, i, j] = (sp.diff(Gamma[l, j, k], coords[i]) -
                                       sp.diff(Gamma[l, i, k], coords[j]) +
                                       sum(Gamma[l, i, m] * Gamma[m, j, k] -
                                           Gamma[l, j, m] * Gamma[m, i, k]
                                           for m in range(n)))

# Step 5: Ricci tensor
Ricci = sp.MutableDenseNDimArray.zeros(n, n)

for i in range(n):
    for j in range(n):
        Ricci[i, j] = sum(Riemann[m, i, m, j] for m in range(n))

# Step 6: Ricci scalar
Ricci_scalar = sum(Ricci[i, j] * metric_inv[i, j] for i in range(n) for j in range(n))

# Step 7: Einstein tensor
Einstein = sp.MutableDenseNDimArray.zeros(n, n)

for i in range(n):
    for j in range(n):
        Einstein[i, j] = Ricci[i, j] - (1 / 2) * metric[i, j] * Ricci_scalar

# Step 8: Pretty-printing Results
print("Metric Tensor (g_ij):")
for i in range(n):
    for j in range(n):
        if metric[i, j] != 0:
            print(f"g[{i}, {j}] =", metric[i, j])

print("\nInverse Metric Tensor (g^ij):")
for i in range(n):
    for j in range(n):
        if metric_inv[i, j] != 0:
            print(f"g[{i}, {j}] =", metric_inv[i, j])

print("\nChristoffel Symbols (Γ^k_ij):")
for k in range(n):
    for i in range(n):
        for j in range(n):
            if Gamma[k, i, j] != 0:
                print(f"Γ[{k}, {i}, {j}] =", Gamma[k, i, j])

print("\nRiemann Tensor (R^l_kij):")
for l in range(n):
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if Riemann[l, k, i, j] != 0:
                    print(f"Riemann[{l}, {k}, {i}, {j}] =", Riemann[l, k, i, j])

print("\nRicci Tensor (R_ij):")
for i in range(n):
    for j in range(n):
        if Ricci[i, j] != 0:
            print(f"R[{i}, {j}] =", Ricci[i, j])

print("\nRicci Scalar (R):", Ricci_scalar)

print("\nEinstein Tensor (G_ij):")
for i in range(n):
    for j in range(n):
        if Einstein[i, j] != 0:
            print(f"G[{i}, {j}] = ", Einstein[i, j])
