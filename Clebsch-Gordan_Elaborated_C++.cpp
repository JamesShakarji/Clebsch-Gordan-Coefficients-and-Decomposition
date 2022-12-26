#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

double clebsch_gordan(double j1, double j2, double j3, double m1, double m2, double m3) {
    // Define the group generators
    double j1, j2, j3;
    // Define the Clebsch–Gordan coefficient
    double CG = sqrt((2 * j3 + 1) / (2 * j1 + 1) / (2 * j2 + 1)) * sum(
        pow(-1, j1 - j2 + k) * sqrt(
            factorial(j1 + j2 - j3 - k) * factorial(j1 - j2 + k) / factorial(j1 + j2 + j3 + 1) / factorial(k) / factorial(j3 - j2 + j1 + k)
        ) * sqrt(
            factorial(j1 + m1) * factorial(j1 - m1) * factorial(j2 + m2) * factorial(j2 - m2)
        ) / sqrt(
            factorial(j3 + m3) * factorial(j3 - m3)
        ) * wigner3j(j1, j2, j3, m1, m2, -m3) * wigner3j(j1, j2, j3, -k, k + m3, m3),
        k, max(0, j3 - j2 + j1), min(j1 + j2 - j3, j1 - m1)
    );
    return CG;
}

int main() {
    // Calculate the Clebsch–Gordan coefficient for j1 = 1, j2 = 1, j3 = 2, m1 = 1, m2 = -1, and m3 = 0
    double CG = clebsch_gordan(1, 1, 2, 1, -1, 0);
    std::cout << CG << std::endl;
    return 0;
}
