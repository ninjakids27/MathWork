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
        double[][] XtX = MatrixOps.Matrix_Mult(MatrixOps.matrixTransposition(xValues), xValues);
        double[][] XtX_inv = MatrixOps.inverseMatrix(XtX);
        double[] XtY = MatrixOps.matrixVectorMult(MatrixOps.matrixTransposition(xValues), yValues);
        double[] coefficients = MatrixOps.matrixVectorMult(XtX_inv, XtY);
        return coefficients;
    }
}
