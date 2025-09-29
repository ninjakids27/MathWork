public class MLOps {
    // makes code more readiable
    private static final double e = Math.E; 
    public static double sigmoid(double x){
        return (1)/(1+Math.pow(e, -x));
    }
    // d stands for derivative
    public static double dSigmoid(double x){
        return (sigmoid(x))*(1-sigmoid(x));
    }
}
