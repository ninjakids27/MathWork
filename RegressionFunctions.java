public class RegressionFunctions {
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
    public static double[] residules(double[][] xValues, double[] yValues, double[] coefficients){
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
     * 
     * @param xValues Data of xValues
     * @param yValues
     * @param timer
     * @return
     */
    public static double[] MLR(double[][] xValues, double[] yValues, boolean timer){
        long startTime = System.nanoTime();
        double[] coefficients = MLR(xValues, yValues);
        long endTime = System.nanoTime();
        long duration = endTime - startTime;
        System.out.println("MLR took: " + String.format("%.9f",MatrixOps.nanoToSeconds(duration)) + " seconds");
        return coefficients;
    }
}
        
