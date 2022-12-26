import java.math.BigDecimal;
import java.math.RoundingMode;

public class ClebschGordan {
    public static BigDecimal clebschGordan(double j1, double j2, double j3, double m1, double m2, double m3) {
        // Define the group generators
        double j1, j2, j3;
        // Define the Clebsch–Gordan coefficient
        BigDecimal CG = sqrt((2 * j3 + 1) / (2 * j1 + 1) / (2 * j2 + 1)).multiply(sum(
            pow(-1, j1 - j2 + k).multiply(sqrt(
                factorial(j1 + j2 - j3 - k).multiply(factorial(j1 - j2 + k))
                    .divide(factorial(j1 + j2 + j3 + 1).multiply(factorial(k)).multiply(factorial(j3 - j2 + j1 + k)))
            )).multiply(sqrt(
                factorial(j1 + m1).multiply(factorial(j1 - m1)).multiply(factorial(j2 + m2)).multiply(factorial(j2 - m2))
            )).divide(sqrt(
                factorial(j3 + m3).multiply(factorial(j3 - m3))
            )).multiply(wigner3j(j1, j2, j3, m1, m2, -m3)).multiply(wigner3j(j1, j2, j3, -k, k + m3, m3)),
            k, max(0, j3 - j2 + j1), min(j1 + j2 - j3, j1 - m1)
        ));
        return CG;
    }

    public static void main(String[] args) {
        // Calculate the Clebsch–Gordan coefficient for j1 = 1, j2 = 1, j3 = 2, m1 = 1, m2 = -1, and m3 = 0
        BigDecimal CG = clebschGordan(1, 1, 2, 1, -1, 0);
        System.out.println(CG);
    }
}
