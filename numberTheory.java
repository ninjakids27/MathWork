public class numberTheory {
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
    private static double rounding(double value, int decimal){
        double scale = Math.pow(10, decimal);
        return Math.round(value * scale) / scale;
    }
    public static double nCr(int n, int r, int decimal){
        return rounding(nCr(n, r), decimal);
    }
    
}
