public class numberTheory {
    /**
     * Rounds a double value to a specified number of decimal places.
     * @param value the value to round
     * @param decimal the number of decimal places
     * @return the rounded value
     */
    private static double rounding(double value, int decimal){
        double scale = Math.pow(10, decimal);
        return Math.round(value * scale) / scale;
    }
    /**
     * Calculates the factorial of a non-negative integer x (x!).
     * @param x the integer to calculate the factorial for
     * @return the factorial of x
     */
    public static long factorial(int x){
        long result = 1;
        for(int i = x; i > 0; i--){
            result *=i;
        }
        return result;
    }
    
    /**
     * Calculates the binomial coefficient "n choose r" (nCr).
     * @param n the total number of items
     * @param r the number of items to choose
     * @return the binomial coefficient
     */
    public static double nCr(int n, int r){
        if (r < 0 || r > n) {
            throw new IllegalArgumentException(ColorText.errorFormat("r must be between 0 and n (inclusive)"));
        }
        // Use a more stable approach to avoid overflow
        double result = 1.0;
        for (int i = 1; i <= r; i++) {
            result *= (n - r + i);
            result /= i;
        }
        return result;
    }
    /**
     * Calculates the binomial coefficient "n choose r" (nCr) and rounds the result to a specified number of decimal places.
     * @param n the total number of items
     * @param r the number of items to choose
     * @param decimal the number of decimal places to round to
     * @return the rounded binomial coefficient
     */
    public static double nCr(int n, int r, int decimal){
        return rounding(nCr(n, r), decimal);
    }
    
}
