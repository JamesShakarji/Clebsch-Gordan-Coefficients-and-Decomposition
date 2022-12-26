import sympy

def clebsch_gordan(j1, j2, j3, m1, m2, m3):
  # Define the group generators
  j1, j2, j3 = sympy.symbols('j1 j2 j3')
  # Define the Clebsch–Gordan coefficient
  CG = sympy.sqrt((2 * j3 + 1) / (2 * j1 + 1) / (2 * j2 + 1)) * sympy.sum(
      (-1) ** (j1 - j2 + k) * sympy.sqrt(
          sympy.factorial(j1 + j2 - j3 - k) * sympy.factorial(j1 - j2 + k) / sympy.factorial(j1 + j2 + j3 + 1) / sympy.factorial(k) / sympy.factorial(j3 - j2 + j1 + k)
      ) * sympy.sqrt(
          sympy.factorial(j1 + m1) * sympy.factorial(j1 - m1) * sympy.factorial(j2 + m2) * sympy.factorial(j2 - m2)
      ) / sympy.sqrt(
          sympy.factorial(j3 + m3) * sympy.factorial(j3 - m3)
      ) * sympy.three_j(j1, j2, j3, m1, m2, -m3) * sympy.three_j(j1, j2, j3, -k, k + m3, m3),
      (k, max(0, j3 - j2 + j1), min(j1 + j2 - j3, j1 - m1))
  )
  return CG

# Calculate the Clebsch–Gordan coefficient for j1 = 1, j2 = 1, j3 = 2, m1 = 1, m2 = -1, and m3 = 0
CG = clebsch_gordan(1, 1, 2, 1, -1, 0)
print(CG)
