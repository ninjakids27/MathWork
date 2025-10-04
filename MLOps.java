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
}
