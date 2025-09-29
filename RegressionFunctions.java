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
    public static void MLR(double[][] xValues, double[] yValues){
        MatrixOps.inverseMatrix(MatrixOps.Matrix_Mult(MatrixOps.matrixTransposition(xValues), xValues))
    }
}
