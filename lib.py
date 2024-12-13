import sympy as sp
import itertools

def get_metric_inv(metric):
    metric_inv = metric.inv()
    return metric_inv

def get_Gamma(metric, metric_inv, coords):
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
    Gamma = sp.simplify(Gamma)
    Gamma = sp.MutableDenseNDimArray(Gamma)
    return Gamma

def get_Riemann(Gamma, coords):
    n = len(coords)
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
    Riemann = sp.simplify(Riemann)
    Riemann = sp.MutableDenseNDimArray(Riemann)
    return Riemann

def get_Ricci(Riemann, coords):
    n = len(coords)
    Ricci = sp.MutableDenseNDimArray.zeros(n, n)

    for i in range(n):
        for j in range(n):
            Ricci[i, j] = sum(Riemann[m, i, m, j] for m in range(n))
    Ricci = sp.simplify(Ricci)
    Ricci = sp.MutableDenseNDimArray(Ricci)
    return Ricci

def get_Ricci_scalar(Ricci, metric_inv, coords):
    n = len(coords)
    Ricci_scalar = sum(Ricci[i, j] * metric_inv[i, j] for i in range(n) for j in range(n))
    Ricci_scalar = sp.simplify(Ricci_scalar)
    return Ricci_scalar

def get_Einstein(Ricci, Ricci_scalar, metric, coords):
    n = len(coords)
    Einstein = sp.MutableDenseNDimArray.zeros(n, n)

    for i in range(n):
        for j in range(n):
            Einstein[i, j] = Ricci[i, j] - (1 / 2) * metric[i, j] * Ricci_scalar
    Einstein = sp.simplify(Einstein)
    Einstein = sp.MutableDenseNDimArray(Einstein)
    return Einstein

def calculate_tensors(metric, coords):
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
    Gamma = sp.simplify(Gamma)
    Gamma = sp.MutableDenseNDimArray(Gamma)

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
    Riemann = sp.simplify(Riemann)
    Riemann = sp.MutableDenseNDimArray(Riemann)

    # Step 4: Ricci tensor
    Ricci = sp.MutableDenseNDimArray.zeros(n, n)

    for i in range(n):
        for j in range(n):
            Ricci[i, j] = sum(Riemann[m, i, m, j] for m in range(n))
    Ricci = sp.simplify(Ricci)
    Ricci = sp.MutableDenseNDimArray(Ricci)

    # Step 5: Ricci scalar
    Ricci_scalar = sum(Ricci[i, j] * metric_inv[i, j] for i in range(n) for j in range(n))
    Ricci_scalar = sp.simplify(Ricci_scalar)

    # Step 6: Einstein tensor
    Einstein = sp.MutableDenseNDimArray.zeros(n, n)

    for i in range(n):
        for j in range(n):
            Einstein[i, j] = Ricci[i, j] - (1 / 2) * metric[i, j] * Ricci_scalar
    Einstein = sp.simplify(Einstein)
    Einstein = sp.MutableDenseNDimArray(Einstein)

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
                print(f"{symbol}{indices} = {value}")

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
                    print(f"{symbol}({i}, {j}) = {value}")

        if is_zero:
            print(f"{symbol} is zero everywhere.")
    
    elif isinstance(tensor, sp.Basic):
        print(f"{symbol} = {tensor}")

    else:
        print(f"Unsupported tensor type: {type(tensor)}")

    print()

def pretty_print(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_scalar, Einstein):
    print_tensor(metric, 'Metric Tensor (g_ij)', 'g')
    print_tensor(metric_inv, 'Inverse Metric Tensor (g^ij)', 'g')
    print_tensor(Gamma, 'Christoffel Symbols (Γ^k_ij)', 'Γ')
    print_tensor(Riemann, 'Riemann Tensor (R^l_kij)', 'R')
    print_tensor(Ricci, 'Ricci Tensor (R_ij)', 'R')
    print_tensor(Ricci_scalar, 'Ricci Scalar (R)', 'R')
    print_tensor(Einstein, 'Einstein Tensor (G_ij)', 'G')

def calculate_print(metric, coords):
    metric_inv, Gamma, Riemann, Ricci, Ricci_curvature, Einstein = calculate_tensors(metric, coords)
    pretty_print(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_curvature, Einstein)

def plt_tensor(tensor, name, symbol):
    latex_string = f'{name}:\n'
    if isinstance(tensor, sp.MutableDenseNDimArray):
        shape = tensor.shape
        is_zero = True

        for indices in itertools.product(*[range(s) for s in shape]):
            value = tensor[indices]
            if value != 0:
                is_zero = False
                latex_string += r"\begin{equation*}" + '\n' + '\t'
                latex_string += f"{symbol}{indices} = {sp.latex(value)}" + '\n'
                latex_string += r"\end{equation*}" + '\n'

        if is_zero:
            latex_string += f"{symbol} is zero everywhere."
    
    elif isinstance(tensor, sp.Matrix):
        is_zero = True
        rows, cols = tensor.shape

        for i in range(rows):
            for j in range(cols):
                value = tensor[i, j]
                if value != 0:
                    is_zero = False
                    latex_string += r"\begin{equation*}" + '\n' + '\t'
                    latex_string += f"{symbol}({i}, {j}) = {sp.latex(value)}" + '\n'
                    latex_string += r"\end{equation*}" + '\n'

        if is_zero:
            latex_string += f"{symbol} is zero everywhere."
    
    elif isinstance(tensor, sp.Basic):
        latex_string += r"\begin{equation*}" + '\n' + '\t'
        latex_string += f"{symbol} = {sp.latex(tensor)}" + '\n'
        latex_string += r"\end{equation*}" + '\n'

    else:
        latex_string += f"Unsupported tensor type: {type(tensor)}"

    latex_string += '\n'
    return latex_string

def pretty_print(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_scalar, Einstein):
    print_tensor(metric, 'Metric Tensor (g_ij)', 'g')
    print_tensor(metric_inv, 'Inverse Metric Tensor (g^ij)', 'g')
    print_tensor(Gamma, 'Christoffel Symbols (Γ^k_ij)', 'Γ')
    print_tensor(Riemann, 'Riemann Tensor (R^l_kij)', 'R')
    print_tensor(Ricci, 'Ricci Tensor (R_ij)', 'R')
    print_tensor(Ricci_scalar, 'Ricci Scalar (R)', 'R')
    print_tensor(Einstein, 'Einstein Tensor (G_ij)', 'G')

def pretty_plt(metric, metric_inv, Gamma, Riemann, Ricci, Ricci_scalar, Einstein):
    latex_string = r"\documentclass{article}" + '\n'
    latex_string += r"\usepackage[fleqn]{amsmath}" + '\n'
    latex_string += r"\usepackage{geometry}" + '\n'
    latex_string += r"\pagenumbering{gobble}" + '\n'
    latex_string += '\t' + r"\geometry{" + '\n'
    latex_string += '\t' + r"left=10mm," + '\n'
    latex_string += '\t' + r"top=10mm," + '\n'
    latex_string += '\t' + r"}" + '\n'
    latex_string += r"\begin{document}" + '\n'

    latex_string += plt_tensor(metric, 'Metric Tensor', 'g')
    latex_string += plt_tensor(metric_inv, 'Inverse Metric Tensor', 'g')
    latex_string += plt_tensor(Gamma, 'Christoffel Symbols', r'\Gamma')
    latex_string += plt_tensor(Riemann, 'Riemann Tensor', 'R')
    latex_string += plt_tensor(Ricci, 'Ricci Tensor', 'R')
    latex_string += plt_tensor(Ricci_scalar, 'Ricci Scalar', 'R')
    latex_string += plt_tensor(Einstein, 'Einstein Tensor', 'G')

    latex_string += r"\end{document}"
    print(latex_string)

    tex_filename = "output.tex"
    with open(tex_filename, 'w') as file:
        file.write(latex_string)
