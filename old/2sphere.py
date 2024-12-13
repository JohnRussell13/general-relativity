import sympy as sp
import itertools

def calculate_tensors(metric):
    # Step 1: Inverese Metric tensor
    metric_inv = metric.inv()

    # Step 2: Christoffel symbols
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

    # Step 3: Riemann tensor
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

    # Step 4: Ricci tensor
    Ricci = sp.MutableDenseNDimArray.zeros(n, n)

    for i in range(n):
        for j in range(n):
            Ricci[i, j] = sum(Riemann[m, i, m, j] for m in range(n))

    # Step 5: Ricci scalar
    Ricci_scalar = sum(Ricci[i, j] * metric_inv[i, j] for i in range(n) for j in range(n))

    # Step 6: Einstein tensor
    Einstein = sp.MutableDenseNDimArray.zeros(n, n)

    for i in range(n):
        for j in range(n):
            Einstein[i, j] = Ricci[i, j] - (1 / 2) * metric[i, j] * Ricci_scalar

    return metric_inv, Gamma, Riemann, Ricci, Ricci_scalar, Einstein

def print_tensor(tensor, name, symbol):
    print(f'{name}:')
    if isinstance(tensor, sp.MutableDenseNDimArray):
        shape = tensor.shape
        is_zero = True

        for indices in itertools.product(*[range(s) for s in shape]):
            value = tensor[indices]
            if value != 0:
                is_zero = False
                print(f"{symbol} {indices}: {value}")

        if is_zero:
            print(f"{symbol} is zero everywhere.")
    
    elif isinstance(tensor, sp.Matrix):
        is_zero = True
        rows, cols = tensor.shape

        for i in range(rows):
            for j in range(cols):
                value = tensor[i, j]
                if value != 0:
                    is_zero = False
                    print(f"{symbol} [{i}, {j}] = {value}")

        if is_zero:
            print(f"{symbol} is zero everywhere.")
    
    elif isinstance(tensor, sp.Basic):
        print(f"{symbol} = {tensor}")

    else:
        print(f"Unsupported tensor type: {type(tensor)}")

    print()

def pretty_print(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_curvature, Einstein):
    print_tensor(metric, 'Metric Tensor (g_ij)', 'g')
    print_tensor(metric_inv, 'Inverse Metric Tensor (g^ij)', 'g')
    print_tensor(Gamma, 'Christoffel Symbols (Γ^k_ij)', 'Γ')
    print_tensor(Riemann, 'Riemann Tensor (R^l_kij)', 'R')
    print_tensor(Ricci, 'Ricci Tensor (R_ij)', 'R')
    print_tensor(Ricci_curvature, 'Ricci Scalar (R)', 'R')
    print_tensor(Einstein, 'Einstein Tensor (G_ij)', 'G')

# Step 0: Needed symbols
R, theta, phi = sp.symbols('R theta phi')
coords = [theta, phi]

# Step 1: Metric tensor
metric = sp.Matrix([
    [R**2, 0],
    [0, R**2 * sp.sin(theta)**2]
])

metric_inv, Gamma, Riemann, Ricci, Ricci_curvature, Einstein = calculate_tensors(metric)

pretty_print(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_curvature, Einstein)

# print("Metric Tensor (g_ij):")
# for i in range(n):
#     for j in range(n):
#         if metric[i, j] != 0:
#             print(f"g[{i}, {j}] =", metric[i, j])
#
# print("\nInverse Metric Tensor (g^ij):")
# for i in range(n):
#     for j in range(n):
#         if metric_inv[i, j] != 0:
#             print(f"g[{i}, {j}] =", metric_inv[i, j])
#
# print("\nChristoffel Symbols (Γ^k_ij):")
# for k in range(n):
#     for i in range(n):
#         for j in range(n):
#             if Gamma[k, i, j] != 0:
#                 print(f"Γ[{k}, {i}, {j}] =", Gamma[k, i, j])
#
# print("\nRiemann Tensor (R^l_kij):")
# for l in range(n):
#     for k in range(n):
#         for i in range(n):
#             for j in range(n):
#                 if Riemann[l, k, i, j] != 0:
#                     print(f"Riemann[{l}, {k}, {i}, {j}] =", Riemann[l, k, i, j])
#
# print("\nRicci Tensor (R_ij):")
# for i in range(n):
#     for j in range(n):
#         if Ricci[i, j] != 0:
#             print(f"R[{i}, {j}] =", Ricci[i, j])
#
# print("\nRicci Scalar (R):", Ricci_scalar)
#
# print("\nEinstein Tensor (G_ij):")
# for i in range(n):
#     for j in range(n):
#         if Einstein[i, j] != 0:
#             print(f"G[{i}, {j}] = ", Einstein[i, j])
