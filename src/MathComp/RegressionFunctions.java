public class RegressionFunctions {
    // tolerance variables for floating point errors 
    private static final double bigEpsilon = 1e11;
    private static final double epsilon = 1e-9;
    /**
     * Checks if two doubles are equal within a defined tolerance (epsilon).
     * @param a first double
     * @param b second double
     * @return true if the numbers are considered equal, false otherwise
     */
    public static boolean isEqual(double a, double b) {
        return Math.abs(a - b) < epsilon;
    }
    /**
     * Uses the defined epsilon for tolerance checking and rounds the input to avoid floating point errors
     * @param input
     * @return
    */
    private static double tolerance(double input){
        double temp = (input*bigEpsilon);
        temp = Math.round(temp);
        return (double)temp/bigEpsilon;
    }
   
    private static void tolernaceArray(double[] array){
        for(int i = 0; i < array.length; i++){
           array[i] = tolerance(array[i]);
        }
    }
    /**
     * Performs standard linear regression on the given x and y values and prints the equation in the form y=mx+b.
     * @param xValues the input x values
     * @param yValues the output y values
     */
    public static void stdLinearRegression(double[] xValues, double[] yValues){
       // prints it out
        double[] coefficients = linearRegression(xValues, yValues);
        System.out.println("y = "+coefficients[0]+"x + "+coefficients[1]);
    }
    
    /**
     * Performs linear regression on the given x and y values using the general MLR routine to ensure consistent handling
     * of the intercept and numerical stability; returns coefficients in order {slope, intercept}.
     * @param xValues the input x values
     * @param yValues the output y values
     * @param onlyLinear if true, only returns the slope (m), if false returns both slope (m) and intercept (b)
     * @return the coefficients of the linear regression (slope and optionally intercept)
     */
    public static double[] linearRegression(double[] xValues, double[] yValues){
        // Build a single-column design matrix (no intercept column) so MLR will add the intercept itself.
        int n = xValues.length;
        double[][] xMatrix = new double[n][1];
        for (int i = 0; i < n; i++) {
            xMatrix[i][0] = xValues[i];
        }
        // MLR returns coefficients in the order [intercept, slope]
        double[] mlrCoeffs = MLR(xMatrix, yValues);
        double intercept = mlrCoeffs.length > 0 ? mlrCoeffs[0] : 0.0;
        double slope = mlrCoeffs.length > 1 ? mlrCoeffs[1] : 0.0;
        return new double[]{slope, intercept};
    }
    public static double[] MLR(double[][] xValues, double[] yValues){
        // Add a column of ones to xValues for the intercept term
        int n = xValues.length;
        int m = xValues[0].length;
        double[][] xValuesWithIntercept = new double[n][m + 1];
        for(int i = 0; i < n; i++){
            xValuesWithIntercept[i][0] = 1; // Intercept term
            for(int j = 0; j < m; j++){
                xValuesWithIntercept[i][j + 1] = xValues[i][j];
            }
        }
        double[][] XtX = MatrixOps.Matrix_Mult(MatrixOps.matrixTransposition(xValuesWithIntercept), xValuesWithIntercept);
        double[][] XtX_inv = MatrixOps.inverseMatrix(XtX);
        double[] XtY = MatrixOps.matrixVectorMult(MatrixOps.matrixTransposition(xValuesWithIntercept), yValues);
        double[] coefficients = MatrixOps.matrixVectorMult(XtX_inv, XtY);
        tolernaceArray(coefficients);
        return coefficients;
    }
    /**
     * Multiple Linear Regression with optional timing for benchmarking.
     * @param xValues Data of xValues
     * @param yValues Data of yValues
     * @param timer Flag to enable timing
     * @return Coefficients of the regression model
     */
    public static double[] MLR(double[][] xValues, double[] yValues, boolean timer){
        long startTime = System.nanoTime();
        double[] coefficients = MLR(xValues, yValues);
        long endTime = System.nanoTime();
        long duration = endTime - startTime;
        System.out.println("MLR took: " + String.format("%.9f",MatrixOps.nanoToSeconds(duration)) + " seconds");
        tolernaceArray(coefficients);
        return coefficients;
    }
    /**
     * Calculates the residuals of a multiple linear regression model.
     * @param xValues the input x values
     * @param yValues the output y values
     * @param coefficients the coefficients of the regression model
     * @return the residuals of the model
     */
    public static double[] residuals(double[][] xValues, double[] yValues, double[] coefficients){
        int n = xValues.length;
        int m = xValues[0].length;
        double[][] xValuesWithIntercept = new double[n][m + 1];
        for(int i = 0; i < n; i++){
            xValuesWithIntercept[i][0] = 1; // Intercept term
            for(int j = 0; j < m; j++){
                xValuesWithIntercept[i][j + 1] = xValues[i][j];
            }
        }
        double[] predictedY = MatrixOps.matrixVectorMult(xValuesWithIntercept, coefficients);
        double[] residuals = new double[n];
        for(int i = 0; i < n; i++){
            residuals[i] = yValues[i] - predictedY[i];
        }
        return residuals;
    }
    /**
     * Performs polynomial regression by transforming the input xValues into polynomial features
     * and then applying multiple linear regression (MLR) to find the coefficients.
     * @param xValues the input x values
     * @param yValues the output y values
     * @param degree the degree of the polynomial
     * @return the coefficients of the polynomial regression
     */
    public static double[] polynomialRegression(double[] xValues, double[] yValues, int degree){
        double[][] xValuesPoly = new double[xValues.length][degree + 1];
        for(int i =0; i < xValues.length; i++){
            for(int j = 0; j <= degree; j++){
                xValuesPoly[i][j] = Math.pow(xValues[i], j);
            }
        }
        double[] coefficients = MLR(xValuesPoly, yValues);
        return coefficients;
    }

    /**
     * Performs power regression by transforming the input x and y values using logarithms
     * and then applying linear regression to find the coefficients.
     * @param x the input x values
     * @param y the output y values
     * @return the coefficients of the power regression
     */
    public static double[] powerRegression(double[] x, double[] y){
        double[] logX = MatrixOps.vectorLog(x);
        double[] logY = MatrixOps.vectorLog(y);
        double[] coefficients = linearRegression(logX, logY);
        coefficients[1] = Math.exp(coefficients[1]);
        return coefficients;
    }

    public static void stdPowerRegression(double[] x, double[] y){
        double[] coefficients = powerRegression(x, y);
        // print it in a format y = ax^b
        System.out.println("y = "+coefficients[1]+"x^"+coefficients[0]);
    }
}
        
