import java.io.BufferedWriter;
import java.io.FileWriter;
public class MLOps {
    
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
     * Build the neural network
     * @param layerSizes array of layer sizes
     * @return 2D array of Neurons representing the network
     * @throws IllegalArgumentException if layerSizes has less than 2 layers
     */
    public static Neuron[][] buildNeuralNetwork(int[] layerSizes) {
        if (layerSizes.length < 2) {
            throw new IllegalArgumentException(ColorText.errorFormat("Neural network must have at least 2 layers."));
        }

        Neuron[][] network = new Neuron[layerSizes.length][];

        for (int i = 0; i < layerSizes.length; i++) {
            network[i] = new Neuron[layerSizes[i]];
            for (int j = 0; j < layerSizes[i]; j++) {
                if (i == 0) {
                    // Input layer neurons
                    network[i][j] = new Neuron(1);
                } else {
                    // Hidden and output layer neurons
                    network[i][j] = new Neuron(layerSizes[i - 1]);
                }
            }
        }
        return network;
    }
    /**
     * Writes neural network structure and weights to a file.
     * first line is structure of the network (e.g., 3,5,2 for a network with 3 input neurons, 5 hidden neurons, and 2 output neurons)
     * or 3,5,5,2 for a network with 3 input neurons, 2 hidden layers of 5 neurons each, and 2 output neurons
     * 
     * this is a save function use the appropiate loading function to read the file back in
     * @param network the neural network to write
     * @param filename the name of the file to write to
     */
    public static void writeNeuralNetworkToFile(Neuron[][] network, String filename) {
        
    }

    /**
     * Loads a neural network from a file.
     * @param filename the name of the file to read from
     * @return the loaded neural network
     */
    public static Neuron[][] loadNeuralNetworkFromFile(String filename) {
        
    }
}
