public class numberTheory {
    /**
     * Calculates the binomial coefficient "n choose r" (nCr).
     * @param n the total number of items
     * @param r the number of items to choose
     * @return the binomial coefficient
     */
    public static double nCr(int n, int r){
        if(r > n) return 0;
        if(r == 0 || r == n) return 1;
        r = Math.min(r, n - r); // Take advantage of symmetry
        long numerator = 1;
        long denominator = 1;
        for(int i = 0; i < r; i++){
            numerator *= (n - i);
            denominator *= (i + 1);
        }
        return (double) numerator / denominator;
    }
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
