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
    public static double residules(double[] function, double[] yValues){
        double[] yPredicted = new double[yValues.length];
        for(int i=0;i<yValues.length;i++){
            yPredicted[i] = function[0]*i + function[1];
        }
        double residules = 0;
        for(int i=0;i<yValues.length;i++){
            residules += Math.pow(yValues[i]-yPredicted[i],2);
        }
        return residules;
    }
}
