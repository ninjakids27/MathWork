
package MLComp;

import Runner.ColorText;
import MLComp.ActivationFunctions_Folder.*;
import MLComp.Optimizers_Folder.*;
import MathComp.Matrix;
import static MathComp.MatrixOps.*;
import java.io.*;

public class MLOps {
    // private static BufferedWriter trainingLogWriter; // Declare trainingLogWriter
    // as a class-level variable
    // tolerance variables for floating point errors
    private static final double bigEpsilon = 1e6;
    private static final double epsilon = 1e-10;

    /**
     * Checks if two doubles are equal within a defined tolerance (epsilon).
     * 
     * @param a first double
     * @param b second double
     * @return true if the numbers are considered equal, false otherwise
     */
    public static boolean isEqual(double a, double b) {
        return Math.abs(a - b) < epsilon;
    }

    /**
     * Uses the defined epsilon for tolerance checking and rounds the input to avoid
     * floating point errors
     * 
     * @param input
     * @return
     */
    private static double tolerance(double input) {
        double temp = (int) Math.round(input * bigEpsilon);
        return (double) temp / bigEpsilon;
    }

    private static void tolernaceArray(double[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = tolerance(array[i]);
        }
    }

    /**
     * Reads a CSV file and returns the values as a double array. in the format as
     * listed in the readme
     * 
     * @param filename
     * @return
     */
    public static int[] readCSV(String filename, int dataNum) {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))){
            int lineNum = 1;
            String line;
            int[] values = new int[785];
            while ((line = br.readLine()) != null) {
                if (lineNum == dataNum) {
                    String[] stringValues = line.split(",");
                    for (int i = 0; i < stringValues.length; i++) {
                        values[i] = Integer.parseInt(stringValues[i]);
                    }
                }
                lineNum++;
            }
            br.close();
            return values;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalAccessError(ColorText.errorFormat("Exception error could not get access to the file: "+filename));
        }
    }

    /**
     * Reads a CSV file and returns the values as a double array. in the format as
     * listed in the readme
     * 
     * @param filename
     * @return
     */
    public static int[][] readCSV(String filename, int start, int end) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            int lineNum = 1;
            String line;
            int[][] values = new int[end - start][785];
            while ((line = br.readLine()) != null) {
                if (lineNum >= start && lineNum < end) {
                    String[] stringValues = line.split(",");
                    values[lineNum - start] = new int[stringValues.length];
                    for (int i = 0; i < stringValues.length; i++) {
                        values[lineNum - start][i] = Integer.parseInt(stringValues[i]);
                    }
                }
                lineNum++;
            }
            br.close();
            return values;
        } catch (Exception e) {
            throw new RuntimeException("IOException occurred while reading CSV file: " + filename, e);
        }
    }

    /**
     * @param network
     * @param input
     * @return
     */
    public static double[] forwardPropagation(Neuron[][] network, double[] input,
            ActivationFunction activationFunction) {
        double[][] inputLayer = new double[network.length][];
        inputLayer[0] = input;
        for (int layer = 1; layer < network.length; layer++) {
            inputLayer[layer] = new double[network[layer].length];
            for (int neuron = 0; neuron < network[layer].length; neuron++) {
                if (layer == network.length - 1) {
                    inputLayer[layer][neuron] = network[layer][neuron].calculateZ(inputLayer[layer - 1]);
                } else {
                    inputLayer[layer][neuron] = network[layer][neuron].activation(inputLayer[layer - 1],
                            activationFunction);
                }
            }
        }
        return softmax(inputLayer[inputLayer.length - 1]);
    }

    /**
     * Forward propagation that returns both pre-activation (z) and post-activation
     * (a) values.
     * 
     * @return double[2][layers][neurons] where [0] = z values, [1] = a values
     *         (activations)
     */
    public static double[][][] forwardPropagationWithZ(Neuron[][] network, double[] input,
            ActivationFunction activationFunction) {
        double[][] zValues = new double[network.length][];
        double[][] aValues = new double[network.length][];

        // Input layer: z = a = input
        zValues[0] = input;
        aValues[0] = input;

        for (int layer = 1; layer < network.length; layer++) {
            zValues[layer] = new double[network[layer].length];
            aValues[layer] = new double[network[layer].length];

            for (int neuron = 0; neuron < network[layer].length; neuron++) {
                // Calculate z = w·a + b
                double z = network[layer][neuron].calculateZ(aValues[layer - 1]);
                zValues[layer][neuron] = z;

                if (layer == network.length - 1) {
                    // Output layer: keep as z (logits), softmax applied later
                    aValues[layer][neuron] = z;
                } else {
                    // Hidden layer: apply activation
                    aValues[layer][neuron] = activationFunction.update(z);
                }
            }
        }

        return new double[][][] { zValues, aValues };
    }

    public static double[][] forwardPropagation(Neuron[][] network, double[] input,
            ActivationFunction activationFunction, boolean returnInputs) {
        double[][] inputLayer = new double[network.length][];
        inputLayer[0] = input;
        for (int layer = 1; layer < network.length; layer++) {
            inputLayer[layer] = new double[network[layer].length];
            for (int neuron = 0; neuron < network[layer].length; neuron++) {
                if (layer == network.length - 1) {
                    inputLayer[layer][neuron] = network[layer][neuron].calculateZ(inputLayer[layer - 1]);
                } else {
                    inputLayer[layer][neuron] = network[layer][neuron].activation(inputLayer[layer - 1],
                            activationFunction);
                }
            }
        }
        return inputLayer;
    }

    public static double[] forwardPropagation(Neuron[][] network, double[] input, boolean timer,
            ActivationFunction activationFunction) {
        double startTime = System.nanoTime();
        double[][] inputLayer = new double[network.length][];
        inputLayer[0] = input;
        for (int layer = 1; layer < network.length; layer++) {
            inputLayer[layer] = new double[network[layer].length];
            for (int neuron = 0; neuron < network[layer].length; neuron++) {
                if (layer == network.length - 1) {
                    inputLayer[layer][neuron] = network[layer][neuron].calculateZ(inputLayer[layer - 1]);
                } else {
                    inputLayer[layer][neuron] = network[layer][neuron].activation(inputLayer[layer - 1],
                            activationFunction);
                }
            }
        }
        double[] result = softmax(inputLayer[inputLayer.length - 1]);
        double endTime = System.nanoTime();
        ColorText.playText(ColorText.dataFormat("Forward Propagation Time: ")
                + ColorText.dataFormat((endTime - startTime) / 1000000.0 + "ms"), 0.1);
        return result;
    }

    /**
     * Mean Squared Error Cost Function
     * 
     * @param output the output from the neural network
     * @param ans    the expected output
     * @return the mean squared error
     */
    public static double lossFunctionMSE(double[] output, double[] ans) {
        double sum = 0.0;
        for (int i = 0; i < output.length; i++) {
            sum += Math.pow((ans[i] - output[i]), 2);
        }
        return sum / output.length;
    }

    /**
     * 
     * 
     * cross entropy loss function
     */
    public static double lossFunctionCSE(double[] output, double[] ans) {
        double sum = 0.0;
        for (int i = 0; i < output.length; i++) {
            sum += -ans[i] * Math.log(output[i] + epsilon);
        }
        return sum;
    }

    /**
     * Build the neural network
     * The structure is that the rows are layers and the columns are neurons in that
     * layer
     * 
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
     * this is a save function use the appropiate loading function to read the file
     * back in
     * 
     * @param network  the neural network to write
     * @param filename the name of the file to write to
     */
    public static void saveNN(Neuron[][] network, String filename) {
        try {
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
     * This uses He uniform initilization which is best suited for ReLU
     * 
     * @param NN
     */
    public static void heinitilizeWeights(Neuron[][] NN) {
        for (int layer = 1; layer < NN.length; layer++) {
            for (int neuron = 0; neuron < NN[layer].length; neuron++) {
                double limit = Math.sqrt(6.0 / NN[layer - 1].length);
                for (int weight = 0; weight < NN[layer][neuron].getWeightsLength(); weight++) {
                    double randomWeight = (Math.random() * 2 * limit) - limit;
                    NN[layer][neuron].setWeight(weight, randomWeight);
                }
            }
        }
    }

    public static void initilizeWeights(Neuron[][] NN) {
        for (int layer = 1; layer < NN.length; layer++) {
            for (int neuron = 0; neuron < NN[layer].length; neuron++) {
                for (int weight = 0; weight < NN[layer][neuron].getWeightsLength(); weight++) {
                    double randomWeight = (Math.random() * 2) - 1;
                    NN[layer][neuron].setWeight(weight, randomWeight);
                }
            }
        }
    }

    /**
     * Loads a neural network from a file.
     * 
     * @param filename the name of the file to read from
     * @return the loaded neural network
     */
    public static Neuron[][] loadNN(String filename) {
        Neuron[][] network = null;
        try {
            FileInputStream file = new FileInputStream(filename);
            ObjectInputStream ois = new ObjectInputStream(file);
            network = (Neuron[][]) ois.readObject();
            ois.close();
            file.close();
            return network;
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
            throw new IllegalArgumentException(ColorText.errorFormat("Model could not be loaded"));
        }
    }

    /**
     * read mnist and forward prop it through
     */

    public static double[] forwardPropagationMNIST(int line, String filename, Neuron[][] NN,
            ActivationFunction activationFunction) {
        int[] input = readCSV(filename, line);
        double[] cleanInput = new double[input.length - 1];
        for (int i = 1; i < input.length; i++) {
            cleanInput[i - 1] = input[i];
        }

        return forwardPropagation(NN, cleanInput, activationFunction);
    }

    public static double[] forwardPropagationMNIST(int line, String filename, Neuron[][] NN,
            ActivationFunction activationFunction, boolean timer) {
        int[] input = readCSV(filename, line);
        double[] cleanInput = new double[input.length - 1];
        for (int i = 1; i < input.length; i++) {
            cleanInput[i - 1] = input[i];
        }

        return forwardPropagation(NN, cleanInput, timer, activationFunction);
    }

    public static double[] softmax(double[] input) {
        // Find max for numerical stability
        double max = input[0];
        for (int i = 1; i < input.length; i++) {
            if (input[i] > max)
                max = input[i];
        }

        // Compute exp(x - max) for numerical stability
        double[] expValues = new double[input.length];
        double sum = 0;
        for (int i = 0; i < input.length; i++) {
            expValues[i] = Math.exp(input[i] - max);
            sum += expValues[i];
        }

        // Normalize
        for (int i = 0; i < expValues.length; i++) {
            expValues[i] /= sum;
        }

        return expValues;
    }

    public static void mnistInterpret(double[] probDistrubution) {
        int index = 0;
        double max = probDistrubution[0];
        for (int i = 0; i < 10; i++) {
            if (max < probDistrubution[i]) {
                index = i;
                max = probDistrubution[i];
            }
        }
        String result = (ColorText.dataFormat("I believe that the number is ") +
                ColorText.returnFormat("" + index) +
                ColorText.dataFormat(" with a confidence of ") +
                ColorText.returnFormat("" + max));
        ColorText.playText(result, 0.05);
    }

    public static void mnistInterpret(double[] probDistrubution, boolean debug) {
        int index = 0;
        double max = probDistrubution[0];
        for (int i = 0; i < 10; i++) {
            if (max < probDistrubution[i]) {
                index = i;
                max = probDistrubution[i];
            }
        }
        String result = (ColorText.dataFormat("I believe that the number is ") +
                ColorText.returnFormat("" + index) +
                ColorText.dataFormat(" with a confidence of ") +
                ColorText.returnFormat("" + max));
        System.out.println(result);
    }

    /**
     * Backpropagation algorithm 
     * @param NN Neural network to train
     * @param filename Dataset to train it on
     * @param learningRate 
     * @param activationFunction 
     * @param activationFunctionDerivative
     * @param OptimizationFunction
     * @param epochAmount The size of the data set to train on so a 60k size dataset and you set it to 10k it will train 10k per epoch of training. 
     */

    public static void backPropagation(Neuron[][] NN, String filename, double learningRate,
            ActivationFunction activationFunction,
            ActivationFunction activationFunctionDerivative,
            Optimizer OptimizationFunction, int epochAmount) {
        int[][] rawInput = readCSV(filename, 0, epochAmount);

        // Normalize inputs to [0, 1] (match test() preprocessing)
        double[][] cleanInputMatrix = new double[rawInput.length][784];
        for (int i = 0; i < rawInput.length; i++) {
            for (int j = 1; j < rawInput[i].length; j++) {
                cleanInputMatrix[i][j - 1] = rawInput[i][j] / 255.0;
            }
        }

        // Iterate through each training example
        for (int trainingRow = 0; trainingRow < cleanInputMatrix.length; trainingRow++) {
            double[] cleanInput = cleanInputMatrix[trainingRow];

            // 2. Forward Pass - get both z (pre-activation) and a (post-activation) values
            double[][][] forwardResult = forwardPropagationWithZ(NN, cleanInput, activationFunction);
            double[][] zValues = forwardResult[0]; // Pre-activation values
            double[][] aValues = forwardResult[1]; // Post-activation values (activations)

            // Get the raw output (logits) from the last layer
            double[] rawOutput = aValues[aValues.length - 1];

            // Apply softmax for prediction and loss calculation
            double[] output = softmax(rawOutput.clone());

            // 3. Create one-hot encoded target
            double[] answer = new double[10];
            answer[rawInput[trainingRow][0]] = 1.0;

            // Backpropagation Start

            // 1. Calculate Deltas
            double[][] deltas = new double[NN.length][];
            for (int layer = NN.length - 1; layer >= 1; layer--) {
                deltas[layer] = new double[NN[layer].length];
                for (int neuron = 0; neuron < NN[layer].length; neuron++) {
                    double delta;
                    if (layer == NN.length - 1) {
                        delta = output[neuron] - answer[neuron];
                    } else {
                        double errorSum = 0.0;
                        Neuron[] nextLayerNeurons = NN[layer + 1];
                        double[] nextLayerDeltas = deltas[layer + 1];
                        for (int k = 0; k < nextLayerNeurons.length; k++) {
                            errorSum += nextLayerDeltas[k] * nextLayerNeurons[k].getWeight(neuron);
                        }
                        delta = errorSum * activationFunctionDerivative.update(zValues[layer][neuron]);
                    }
                    deltas[layer][neuron] = delta;
                }
            }
            // calculate gradients and update weights using the optimizer
            for (int layer = NN.length - 1; layer >= 1; layer--) {
                for (int neuron = 0; neuron < NN[layer].length; neuron++) {
                    double delta = deltas[layer][neuron];
                    double[] inputs = aValues[layer - 1];
                    double[] gradient = new double[inputs.length + 1];

                    for (int w = 0; w < inputs.length; w++) {
                        gradient[w] = inputs[w] * delta;
                    }
                    gradient[gradient.length - 1] = delta;

                    OptimizationFunction.update(NN[layer][neuron], gradient, learningRate);
                }
            }
            // Progress reporting every 1000 examples
            if ((trainingRow + 1) % epochAmount == 0) {
                // Compute accuracy on a fixed-size recent subset to avoid O(n^2) logging cost
                int windowSize = 1000;
                int startIdx = Math.max(0, (trainingRow + 1) - windowSize);
                double correctPredictions = 0;
                for (int i = startIdx; i <= trainingRow; i++) {
                    double[] predOutput = forwardPropagation(NN, cleanInputMatrix[i], activationFunction);
                    int predictedLabel = (int) argmax(predOutput);
                    if (predictedLabel == rawInput[i][0]) {
                        correctPredictions++;
                    }
                }
                int totalEvaluated = trainingRow - startIdx + 1;
                double accuracy = (correctPredictions / totalEvaluated) * 100.0;
                String returnString = ColorText.dataFormat("Training Progress: ") +
                        ColorText.returnFormat((trainingRow + 1) + "/" + cleanInputMatrix.length) +
                        ColorText.dataFormat(" - Accuracy: ") +
                        ColorText.returnFormat(String.format("%.2f", accuracy) + "%") +
                        ColorText.dataFormat(" - Current Loss: ") +
                        ColorText.returnFormat(String.format("%.6f", lossFunctionCSE(output, answer)));
                System.out.println(returnString);
            }
        }
    }

    public static void training(Neuron[][] NN, String filename, double learningRate, int epochs,
            ActivationFunction activationFunction, ActivationFunction activationFunctionDerivative,
            Optimizer OptimizationFunction, int epochAmount) {
        for (int epoch = 0; epoch < epochs; epoch++) {
            System.out
                    .println(ColorText.dataFormat("Epoch ") + ColorText.returnFormat("" + (epoch + 1) + "/" + epochs));
            backPropagation(NN, filename, learningRate, activationFunction, activationFunctionDerivative,
                    OptimizationFunction, epochAmount);
        }
    }

    public static void training(Neuron[][] NN, String filename, double learningRate, int epochs,
            ActivationFunction activationFunction, ActivationFunction activationFunctionDerivative,
            Optimizer OptimizationFunction, int epochAmount, boolean test) {
        for (int epoch = 0; epoch < epochs; epoch++) {
            System.out
                    .println(ColorText.dataFormat("Epoch ") + ColorText.returnFormat("" + (epoch + 1) + "/" + epochs));
            backPropagation(NN, filename, learningRate, activationFunction, activationFunctionDerivative,
                    OptimizationFunction, epochAmount);
            if (test) {
                double accuracy = (double) test(NN, "MNIST_CSV/mnist_test.csv", activationFunction,
                        activationFunctionDerivative) / 10000.0;
                System.out.println(ColorText.dataFormat("Current Accuracy: ")
                        + ColorText.returnFormat(String.format("%.2f", accuracy * 100.0) + "%"));
            }
        }
    }

    public static void training(Neuron[][] NN, String filename, double learningRate, double accuracy,
            ActivationFunction activationFunction, ActivationFunction activationFunctionDerivative,
            Optimizer OptimizationFunction, int epochAmount) {
        double currentAccuracy = 0.0;
        int epoch = 0;
        while (currentAccuracy < accuracy) {
            
            System.out.println(ColorText.dataFormat("Epoch ") + ColorText.returnFormat("" + (epoch + 1)));
            backPropagation(NN, filename, learningRate, activationFunction, activationFunctionDerivative,
                    OptimizationFunction, epochAmount);
            currentAccuracy = (double) test(NN, "MNIST_CSV/mnist_test.csv", activationFunction,
                    activationFunctionDerivative) / 10000.0;
            System.out.println(ColorText.dataFormat("Current Accuracy: ")
                    + ColorText.returnFormat(String.format("%.2f", currentAccuracy * 100.0) + "%"));
            epoch++;
        }
    }

    public static void training(Neuron[][] NN, String filename, double learningRate, double accuracy,
            ActivationFunction activationFunction, ActivationFunction activationFunctionDerivative,
            Optimizer OptimizationFunction, int epochAmount, boolean noTest) {
        double currentAccuracy = 0.0;
        int epoch = 0;
        while (currentAccuracy < accuracy) {
            System.out.println(ColorText.dataFormat("Epoch ") + ColorText.returnFormat("" + (epoch + 1)));
            backPropagation(NN, filename, learningRate, activationFunction, activationFunctionDerivative,
                    OptimizationFunction, epochAmount);

            epoch++;
        }
    }

    public static int test(Neuron[][] NN, String filename, ActivationFunction activationFunction,
            ActivationFunction activationFunctionDerivative) {
        int[][] rawInput = readCSV(filename, 0, 10_000); // Test on 10,000 examples!
        // Normalize inputs
        double[][] cleanInputMatrix = new double[rawInput.length][784];
        for (int i = 0; i < rawInput.length; i++) {
            for (int j = 1; j < rawInput[i].length; j++) {
                cleanInputMatrix[i][j - 1] = rawInput[i][j] / 255.0;
            }
        }
        // do forward pass and count the correct outputs
        int correct = 0;
        for (int i = 0; i < rawInput.length; i++) {
            double[] cleanInput = new double[784];
            for (int j = 1; j < rawInput[i].length; j++) {
                cleanInput[j - 1] = rawInput[i][j] / 255.0;
            }
            double[] output = forwardPropagation(NN, cleanInput, activationFunction);
            int predictedLabel = (int) argmax(output);
            if (predictedLabel == rawInput[i][0]) {
                correct++;
            }
        }
        return correct;
    }
}
