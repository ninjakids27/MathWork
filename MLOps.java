import java.io.*;
public class MLOps {
    /**
     * Reads a CSV file and returns the values as a double array. in the format as listed in the readme
     * @param filename
     * @return
     */
    public static double[] readCSV(String filename, int bob){
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            int lineNum = 1;
            String line;
            double[] values = new double[785];
            while((line = br.readLine()) != null){
                if(lineNum == bob){
                    String[] stringValues = line.split(",");
                    for (int i = 0; i < stringValues.length; i++) {
                        values[i] = Double.parseDouble(stringValues[i]);
                    }
                }
                lineNum++;
            }
            br.close();
            return values;
        } catch (Exception e) {
            System.out.println("IOException return is garbage");
            return new double[0];
        }
    }

    /**
     * 
     * @param network
     * @param input
     * @return
     */
    public static double[] forwardPropagation(Neuron[][] network, double[] input) {
        double[][] inputLayer = new double[network.length][];
        inputLayer[0] = input;
        for(int layer = 1; layer < network.length; layer++) {
            inputLayer[layer] = new double[network[layer].length];
            for(int neuron = 0; neuron < network[layer].length; neuron++){
                inputLayer[layer][neuron] = network[layer][neuron].activation(inputLayer[layer-1]);    
            }
        }
        return inputLayer[inputLayer.length - 1];
    }

    /**
     * Mean Squared Error Cost Function
     * @param output the output from the neural network
     * @param ans the expected output
     * @return the mean squared error
     */
    public static double costFunction(double[] output, double[] ans){
        if(output.length != ans.length){
            throw new IllegalArgumentException(ColorText.errorFormat("Output and answer arrays must be of the same length."));
        }
        double sum = 0.0;
        for(int i = 0; i < output.length; i++){
            sum += Math.pow(output[i] - ans[i], 2);
        }
        return sum / output.length;
    }
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
     * The structure is that the rows are layers and the columns are neurons in that layer 
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
     * Serilizes NN to a file
     * 
     * this is a save function use the appropiate loading function to read the file back in
     * @param network the neural network to write
     * @param filename the name of the file to write to
     */
    public static void saveNN(Neuron[][] network, String filename) {
        try{
            FileOutputStream file = new FileOutputStream(filename);
            ObjectOutputStream oos = new ObjectOutputStream(file);
            oos.writeObject(network);
            oos.close();
            file.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * This uses He initilization which is best suited for ReLU
     * @param NN
     */
    public static void initilizeWeights(Neuron[][] NN){

    }
    /**
     * Loads a neural network from a file.
     * @param filename the name of the file to read from
     * @return the loaded neural network
     */
    public static Neuron[][] loadNN(String filename) {
        try{
            FileInputStream file = new FileInputStream(filename);
            ObjectInputStream ois = new ObjectInputStream(file);
            Neuron[][] network = null;
            network = (Neuron[][]) ois.readObject();
            ois.close();
            file.close();
            return network;
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
            return null;
        }
    }
}
