public class sandbox {
    public static void main(String[] args){
        double[] x = {1.227,2.12,6.283,4.24,3.347};
        double[] y = {0.2,0.4,0.8,0.6,0.5};
        
        // MatrixOps.vectorPow(x, 0.0001);

        RegressionFunctions.stdLinearRegression(x, y);
    }
}
