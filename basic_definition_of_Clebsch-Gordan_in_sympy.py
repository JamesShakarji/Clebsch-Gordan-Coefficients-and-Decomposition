import sympy

# Define the symbols for the angular momenta
j1, j2, m1, m2 = sympy.symbols('j1 j2 m1 m2')

# Define the Clebsch–Gordan coefficients
cg = sympy.clebsch_gordan(j1, j2, m1, m2)

# Print the result
print(cg)
