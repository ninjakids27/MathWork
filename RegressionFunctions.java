public class RegressionFunctions {
    /**
     * Performs standard linear regression on the given x and y values and prints the equation in the form y=mx+b.
     * @param xValues the input x values
     * @param yValues the output y values
     */
    public static void stdLinearRegression(double[] xValues, double[] yValues){
       int n = xValues.length;
        double[] xValuesCopy = xValues.clone();
        double m = (n*StatOps.sum(StatOps.vectorProduct(xValues, yValues))-StatOps.sum(xValues)*StatOps.sum(yValues))
                    /
                    (n*StatOps.sum(StatOps.vectorPow(xValuesCopy, 2))-Math.pow((StatOps.sum(xValues)),2));
                    
                    
        double b = (StatOps.sum(yValues)-m*StatOps.sum(xValues))/n;
        System.out.println("y="+m+"x+"+b);
    }
    
    public static double[] linearRegression(double[] xValues, double[] yValues,boolean onlyLinear){
       int n = xValues.length;
        double[] xValuesCopy = xValues.clone();
        double m = (n*StatOps.sum(StatOps.vectorProduct(xValues, yValues))-StatOps.sum(xValues)*StatOps.sum(yValues))
                    /
                    (n*StatOps.sum(StatOps.vectorPow(xValuesCopy, 2))-Math.pow((StatOps.sum(xValues)),2));
                    
        if(!onlyLinear){
            double b = (StatOps.sum(yValues)-m*StatOps.sum(xValues))/n;
            double[] a = {m,b};
            return a;
        }

        double[] a = {m};
        return a;
    }
    /**
     * Performs linear regression on the given x and y values. It uses the formula instead of actual linear algebra
     * @param xValues the input x values
     * @param yValues the output y values
     * @param onlyLinear if true, only returns the slope (m), if false returns both slope (m) and intercept (b) (returns array of size 1 or 2)
     * @return the coefficients of the linear regression (slope and intercept)
     */
    public static double[] linearRegression(double[] xValues, double[] yValues){
       int n = xValues.length;
        double[] xValuesCopy = xValues.clone();
        double m = (n*StatOps.sum(StatOps.vectorProduct(xValues, yValues))-StatOps.sum(xValues)*StatOps.sum(yValues))
                    /
                    (n*StatOps.sum(StatOps.vectorPow(xValuesCopy, 2))-Math.pow((StatOps.sum(xValues)),2));
                    
                    
        double b = (StatOps.sum(yValues)-m*StatOps.sum(xValues))/n;
        
        double[] a = {m,b};
        return a;
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
        double[] logX = StatOps.vectorLog(x);
        double[] logY = StatOps.vectorLog(y);
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
        
