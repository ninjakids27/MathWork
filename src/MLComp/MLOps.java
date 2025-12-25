import java.io.*;
public class MLOps {
    // tolerance variables for floating point errors 
    private static final double bigEpsilon = 1e6;
    private static final double epsilon = 1e-10;
    /**
     * Checks if two doubles are equal within a defined tolerance (epsilon).
     * @param a first double
     * @param b second double
     * @return true if the numbers are considered equal, false otherwise
     */
    public static boolean isEqual(double a, double b) {
        return Math.abs(a - b) < epsilon;
    }
    /**
     * Uses the defined epsilon for tolerance checking and rounds the input to avoid floating point errors
     * @param input
     * @return
    */
    private static double tolerance(double input){
        double temp = (int) Math.round(input*bigEpsilon);
        return (double)temp/bigEpsilon;
    }
   
    private static void tolernaceArray(double[] array){
        for(int i = 0; i < array.length; i++){
           array[i] = tolerance(array[i]);
        }
    }
    
    /**
     * Reads a CSV file and returns the values as a double array. in the format as listed in the readme
     * @param filename
     * @return
     */
    public static int[] readCSV(String filename, int dataNum){
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            int lineNum = 1;
            String line;
            int[] values = new int[785];
            while((line = br.readLine()) != null){
                if(lineNum == dataNum){
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
            System.out.println("IOException return is garbage");
            return new int[0];
        }
    }

    /**
     *  
     * @param network
     * @param input
     * @return
     */
    public static double[] forwardPropagation(Neuron[][] network, double[] input, ActivationFunction activationFunction) {
        double[][] inputLayer = new double[network.length][];
        inputLayer[0] = input;
        for(int layer = 1; layer < network.length; layer++) {
            inputLayer[layer] = new double[network[layer].length];
            for(int neuron = 0; neuron < network[layer].length; neuron++){
                inputLayer[layer][neuron] = network[layer][neuron].activation(inputLayer[layer-1],activationFunction);    
            }
        }
        return softmax(inputLayer[inputLayer.length - 1]);
    }

    public static double[][] forwardPropagation(Neuron[][] network, double[] input, ActivationFunction activationFunction,boolean returnInputs) {
        double[][] inputLayer = new double[network.length][];
        inputLayer[0] = input;
        for(int layer = 1; layer < network.length; layer++) {
            inputLayer[layer] = new double[network[layer].length];
            for(int neuron = 0; neuron < network[layer].length; neuron++){
                inputLayer[layer][neuron] = network[layer][neuron].activation(inputLayer[layer-1],activationFunction);    
            }
        }
        return inputLayer;
    }

    
    public static double[] forwardPropagation(Neuron[][] network, double[] input,boolean timer, ActivationFunction activationFunction) {
        double startTime = System.nanoTime();
        double[][] inputLayer = new double[network.length][];
        inputLayer[0] = input;
        for(int layer = 1; layer < network.length; layer++) {
            inputLayer[layer] = new double[network[layer].length];
            for(int neuron = 0; neuron < network[layer].length; neuron++){
                inputLayer[layer][neuron] = network[layer][neuron].activation(inputLayer[layer-1],activationFunction);    
            }
        }
        double[] result = softmax(inputLayer[inputLayer.length - 1]);
        double endTime = System.nanoTime();
        ColorText.playText(ColorText.dataFormat("Forward Propagation Time: ")+ColorText.dataFormat((endTime-startTime)/1000000.0+"ms"),0.1);
        return result;
    }
    
    /**
     * Mean Squared Error Cost Function
     * @param output the output from the neural network
     * @param ans the expected output
     * @return the mean squared error
     */
    public static double lossFunctionMSE(double[] output, double[] ans){
        double sum = 0.0;
        for(int i = 0; i < output.length; i++){
            sum += Math.pow((ans[i]-output[i]), 2);
        }
        return sum/output.length;
    }

    /**
     * 
     * 
     * cross entropy loss function
     */
    public static double lossFunctionCSE(double[] output, double[] ans){
        double sum = 0.0;
        for(int i = 0; i < output.length; i++){
            sum += -ans[i]*Math.log(output[i]+epsilon);
        }
        return sum;
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
     * This uses He uniform initilization which is best suited for ReLU
     * @param NN
     */
    public static void heinitilizeWeights(Neuron[][] NN){
        for(int layer = 1; layer < NN.length; layer++){
            for(int neuron = 0; neuron < NN[layer].length; neuron++){
                double limit = Math.sqrt(6.0/NN[layer-1].length);
                for(int weight = 0; weight < NN[layer][neuron].getWeightsLength(); weight++){
                    double randomWeight = (Math.random()*2*limit)-limit;
                    NN[layer][neuron].setWeight(weight, randomWeight);
                }
            }
        }
    }

    public static void initilizeWeights(Neuron[][] NN){
        for(int layer = 1; layer < NN.length; layer++){
            for(int neuron = 0; neuron < NN[layer].length; neuron++){
                for(int weight = 0; weight < NN[layer][neuron].getWeightsLength(); weight++){
                    double randomWeight = (Math.random()*2)-1;
                    NN[layer][neuron].setWeight(weight, randomWeight);
                }
            }
        }
    }
    /**
     * Loads a neural network from a file.
     * @param filename the name of the file to read from
     * @return the loaded neural network
     */
    public static Neuron[][] loadNN(String filename) {
        try{
            Neuron[][] network = null;
            FileInputStream file = new FileInputStream(filename);
            ObjectInputStream ois = new ObjectInputStream(file);
            network = (Neuron[][]) ois.readObject();
            ois.close();
            file.close();
            return network;
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
            return null;
        }
    }

    /**
     * read mnist and forward prop it through
     */

    public static double[] forwardPropagationMNIST(int line,String filename, Neuron[][] NN,ActivationFunction activationFunction){
        int[] input = readCSV(filename, line);
        double[] cleanInput = new double[input.length-1];
        for(int i = 1; i < input.length; i++){
            cleanInput[i-1] = input[i];
        }

        return forwardPropagation(NN, cleanInput,activationFunction);
    }

    public static double[] forwardPropagationMNIST(int line,String filename, Neuron[][] NN,ActivationFunction activationFunction,boolean timer){
        int[] input = readCSV(filename, line);
        double[] cleanInput = new double[input.length-1];
        for(int i = 1; i < input.length; i++){
            cleanInput[i-1] = input[i];
        }

        return forwardPropagation(NN, cleanInput,timer,activationFunction);
    }

    public static double[] softmax(double[] input){
        MatrixOps.vectorPow(input, Math.E);
        double sum = 0;
        for(double num : input){
            sum+=num;
        }
        MatrixOps.vectorMult(input, 1.0/sum);
        tolernaceArray(input);
        return input;
    }

    public static void mnistInterpret(double[] probDistrubution){
        int index = 0;
        double max = probDistrubution[0]; 
        for(int i = 0; i < 10; i++){
            if(max < probDistrubution[i]){
                index = i;
                max = probDistrubution[i];
        }
    }
    String result = (
        ColorText.dataFormat("I believe that the number is ")+
        ColorText.returnFormat(""+index)+
        ColorText.dataFormat(" with a confidence of ")+
        ColorText.returnFormat(""+max)
    );
    ColorText.playText(result,0.05);
    }

    public static void mnistInterpret(double[] probDistrubution, boolean debug){
        int index = 0;
        double max = probDistrubution[0]; 
        for(int i = 0; i < 10; i++){
            if(max < probDistrubution[i]){
                index = i;
                max = probDistrubution[i];
        }
    }
    String result = (
        ColorText.dataFormat("I believe that the number is ")+
        ColorText.returnFormat(""+index)+
        ColorText.dataFormat(" with a confidence of ")+
        ColorText.returnFormat(""+max)
    );
    System.out.println(result);
    }

    /**
     * 
     * @param NN
     * @param line
     * @param filename
     * @param learningRate
     * @param activationFunction
     * @param activationFunctionDerivative
     */
    public static void backPropagation(Neuron[][] NN, int line, String filename, double learningRate, ActivationFunction activationFunction, ActivationFunction activationFunctionDerivative, Optimizer OptimizationFunction) {
    // 1. Data Prep & Normalization
    int[] rawInput = readCSV(filename, line);
    double[] answer = new double[10];
    answer[rawInput[0]] = 1.0; // One-hot encoding

    // Normalize inputs (crucial for DNNs to prevent saturation)
    double[] cleanInput = new double[rawInput.length - 1];
    for (int i = 1; i < rawInput.length; i++) {
        cleanInput[i - 1] = rawInput[i] / 255.0; 
    }

    // 2. Forward Pass (Capture all intermediate activations)
    double[] output = forwardPropagation(NN, cleanInput, activationFunction);
    double[][] layerActivations = forwardPropagation(NN, cleanInput, activationFunction, true);

    // 3. Output Layer Error
    // error = (output - target) * deriv(output)
    double[] currentError = new double[output.length];
    for (int i = 0; i < output.length; i++) {
        currentError[i] = (output[i] - answer[i]) * activationFunctionDerivative.update(output[i]);
    }

    // 4. Backpropagate Errors through Hidden Layers
    // Store errors for every layer so we can use them for weight updates
    double[][] allLayerErrors = new double[NN.length][];
    allLayerErrors[NN.length - 1] = currentError;

    for (int layer = NN.length - 2; layer > 0; layer--) {
        double[] nextLayerError = allLayerErrors[layer + 1];
        double[] thisLayerError = new double[NN[layer].length];

        for (int i = 0; i < NN[layer].length; i++) {
            double errorSum = 0.0;
            for (int j = 0; j < NN[layer + 1].length; j++) {
                errorSum += NN[layer + 1][j].getWeight(i) * nextLayerError[j];
            }
            thisLayerError[i] = errorSum * activationFunctionDerivative.update(layerActivations[layer][i]);
        }
        allLayerErrors[layer] = thisLayerError;
    }

    // 5. Update Weights using Gradients
    // Gradient = previous_layer_activation * current_layer_error
    for (int layer = 1; layer < NN.length; layer++) {
        double[] layerError = allLayerErrors[layer];
        double[] prevLayerActivations = (layer == 1) ? cleanInput : layerActivations[layer - 1];

        for (int neuron = 0; neuron < NN[layer].length; neuron++) {
            double[] gradient = new double[NN[layer][neuron].getWeightsLength()];
            for (int w = 0; w < gradient.length; w++) {
                // Weight update formula
                gradient[w] = prevLayerActivations[w] * layerError[neuron] * learningRate;
            }
            OptimizationFunction.update(NN[layer][neuron], gradient);
        }
    }

    // Progress Reporting
    if (true) { // showProgress
        double loss = lossFunctionMSE(output, answer);
        int predicted = 0;
        double max = output[0];
        for (int i = 1; i < output.length; i++) {
            if (output[i] > max) {
                max = output[i];
                predicted = i;
            }
        }
    System.out.printf("%-28s | %-15s | %-12s | %-12s | %-12s%n",
        ColorText.dataFormat("Data Line: ") + ColorText.returnFormat(""+line),
        ColorText.dataFormat("Correctness: ") + ColorText.returnFormat(""+(predicted == rawInput[0])),
        ColorText.dataFormat("Loss: ") + ColorText.returnFormat(String.format("%.6f", loss)),
        ColorText.dataFormat("Predicted: ") + ColorText.returnFormat(""+predicted),
        ColorText.dataFormat("Actual: ") + ColorText.returnFormat(""+rawInput[0])
    );
    }
}
    
    public static void training(Neuron[][] NN, String filename, double learningRate, int epochs, ActivationFunction activationFunction, ActivationFunction activationFunctionDerivative, Optimizer OptimizationFunction){
        for(int epoch = 0; epoch < epochs; epoch++){
            ColorText.playText(ColorText.dataFormat("Epoch ")+ColorText.returnFormat(""+(epoch+1)+"/"+epochs),0.1);

            for(int line = 1; line < 60000; line++){
                backPropagation(NN, line, filename, learningRate, activationFunction, activationFunctionDerivative, OptimizationFunction);
            }
        }
    }
}
