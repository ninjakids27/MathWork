package MLComp.ActivationFunctions_Folder;

public class ActivationFunctions {
    /**
     * ReLU activation function (see notes for specifics)
     * @param x input value
     * @return output value
     */
    public static double reLU(double x) {
        return Math.max(0, x);
    }
    /**
     * Derivative of ReLU activation function
     * @param x input value
     * @return derivative value
     */
    public static double reLUDerivative(double x) {
        return x > 0 ? 1 : 0;
    }
    /**
     * Sigmoid activation function
     * @param x input value
     * @return output value
     */
    public static double sigmoid(double x) {
        return 1 / (1 + Math.exp(-x));
    }
    /**
     * Derivative of Sigmoid activation function
     * @param x input value
     * @return derivative value
     */
    public static double sigmoidDerivative(double x) {
        double sig = sigmoid(x);
        return sig * (1 - sig);
    }

    /**
     * leaky relu activation function
     * @param x input value
     */
    public static double leakyReLU(double x) {
        return x > 0 ? x : 0.01 * x;
    }
    /**
     * Derivative of leaky relu activation function
     * @param x input value
     */
    public static double leakyReLUDerivative(double x) {
        return x > 0 ? 1 : 0.01;
    }

    public static double tanh(double x) {
        return Math.tanh(x);
    }

    public static double tanhDerivative(double x) {
        double t = tanh(x);
        return 1 - t * t;
    }
}
