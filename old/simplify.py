import sympy as sp

# Declare the variables
m, r, theta = sp.symbols('m r theta')

# Define the full expression for G(0,0)
expr = -1.0*m**2/r**4 - 1.0*m**2*(-2*m + r)**2/(r**6*(-2*m/r + 1)**2) + 1.0*m/r**3 - 1.0*m*(-2*m + r)/r**4 - (1.0*m/r - 0.5)*(-r*(-1.0*m**2/r**4 - 1.0*m**2*(-2*m + r)**2/(r**6*(-2*m/r + 1)**2) + 1.0*m/r**3 - 1.0*m*(-2*m + r)/r**4)/(-2*m + r) + (-2*m + r)*(-1.0*m**2/(r**2*(-2*m + r)**2) - 1.0*m**2/(r**4*(-2*m/r + 1)**2) + 1.0*m/(r*(-2*m + r)**2) + 1.0*m/(r**2*(-2*m + r)) - 2.0*m*(-2*m + r)/(r**4*(-2*m/r + 1)**2))/r + (-1.0*m*sp.sin(theta)**2/r + 1.0*m*(-2*m + r)**2*sp.sin(theta)**2/(r**3*(-2*m/r + 1)**2))/(r**2*sp.sin(theta)**2) + (1.0*m*(2.0*m - 1.0*r)/(r*(-2*m + r)) - 1.0*m*(-2*m + r)*(2.0*m - 1.0*r)/(r**3*(-2*m/r + 1)**2))/r**2)

# Simplify the expression
simplified_expr = sp.simplify(expr)

# Optionally, expand and then factor to get even more simplifications
expanded_expr = sp.expand(simplified_expr)
factored_expr = sp.factor(expanded_expr)

# Print the simplified, expanded, and factored expressions
print("Simplified:", simplified_expr)
print("Expanded:", expanded_expr)
print("Factored:", factored_expr)
