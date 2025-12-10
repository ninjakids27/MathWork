/**
 * This class represents a single neuron in a neural network.
 * It contains weights and a bias, and provides methods to get and set these values.
 */

public class Neuron implements java.io.Serializable {
    private double[] weights;
    private double bias;
    /**
     * Constructor to initialize the neuron with a given number of inputs.
     * @param inputSize
     */
    public Neuron(int inputSize) {
        weights = new double[inputSize];
        bias = 0.0;
    }

    /**
     * Constructor to initialize the neuron with given weights and bias.
     * @param weights the weights array
     * @param bias the bias value
     */
    public Neuron(double[] weights, double bias) {
        this.weights = weights;
        this.bias = bias;
    }

    /**
     * Gets the weights of the neuron.
     * @return the weights array
     */
    public double[] getWeights() {
        return weights;
    }
    /**
     * Gets the bias of the neuron.
     * @return the bias value
     */
    public double getBias() {
        return bias;
    }
    public int getWeightsLength(){
        return this.weights.length;
    }
    // methods to update the weights and bias
    /**
     * Sets the weights of the neuron.
     * @param weights the new weights array
     * @throws IllegalArgumentException if the length of the weights array does not match the input size
     */
    public void setWeights(double[] weights) {
        if(weights.length != this.weights.length) {
            throw new IllegalArgumentException("Weights array length must match input size.");
        }
        this.weights = weights;
    }


    public void setWeight(int index, double weight){
        try {
            // System.out.println(weight);
            weights[index] = weight;
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new IllegalArgumentException(ColorText.errorFormat("Cannot put the weight in the weights for the neuron"));
        }
    }
    /**
     * Sets the bias of the neuron.
     * @param bias the new bias value
     */
    public void setBias(double bias) {
        this.bias = bias;
    }
    /**
     * Prints the weights and bias of the neuron to the console.
     */
    public void print() {
        System.out.print("Weights: ");
        for(double weight : weights) {
            System.out.print(weight + " ");
        }
        System.out.println("\nBias: " + bias);
    }
    /**
     * Computes the activation of the neuron given an array of inputs. using ReLU activation function.
     * @param inputs
     * @return
     */
    public double activation(double[] inputs){
        return MLOps.sigmoid(MatrixOps.dot_product(weights, inputs)+bias);
    }

    public double softmax(){
        
    }
}
