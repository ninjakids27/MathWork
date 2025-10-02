public class Neuron {
    private double[] weights;
    private double bias;
    
    public Neuron(int inputSize) {
        weights = new double[inputSize];
        bias = 0.0;
    }

    public double[] getWeights() {
        return weights;
    }

    public double getBias() {
        return bias;
    }

    // methods to update the weights and bias
    public void setWeights(double[] weights) {

        if(weights.length != this.weights.length) {
            throw new IllegalArgumentException("Weights array length must match input size.");
        }
        this.weights = weights;
    }


    public void setBias(double bias) {
        this.bias = bias;
    }
}
