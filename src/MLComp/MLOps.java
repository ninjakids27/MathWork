import java.io.*;
import java.util.function.Function;
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
    public static double[] forwardPropagation(Neuron[][] network, double[] input, Function<Double, Double> activationFunction) {
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
    public static double[] forwardPropagation(Neuron[][] network, double[] input,boolean timer, Function<Double, Double> activationFunction) {
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
    public static double lossFunction(double[] output, double[] ans){
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
     * This uses He uniform initilization which is best suited for ReLU
     * @param NN
     */
    public static void initilizeWeights(Neuron[][] NN){
        // go through each hidden layer
        for(int i = 1; i < NN.length; i++){
            Neuron[] layer = NN[i];
            int n_low   = NN[i-1].length; 
            double low  = Math.sqrt(6.0/n_low);
            // initilize the weights based on the neurons per layer
            for(Neuron neuron : layer){
                // iterate through the weights sampling from a uniform distribution
                for(int weight = 0; weight < neuron.getWeightsLength(); weight++){
                    // he uniform is U(-L,L)
                     neuron.setWeight(weight, -1*low+Math.random()*(2*low));
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

    /**
     * read mnist and forward prop it through
     */

    public static double[] forwardPropagationMNIST(int line,String filename, Neuron[][] NN,Function<Double, Double> activationFunction){
        int[] input = readCSV(filename, line);
        double[] cleanInput = new double[input.length-1];
        for(int i = 1; i < input.length; i++){
            cleanInput[i-1] = input[i];
        }

        return forwardPropagation(NN, cleanInput,activationFunction );
    }

    public static double[] forwardPropagationMNIST(int line,String filename, Neuron[][] NN,Function<Double, Double> activationFunction,boolean timer){
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
     */
    public static void backPropagation(Neuron[][] NN){
        

    }
    /**
     * 
     */
    public static void ADAM(){

    }



    public static void SGD(){

    }

    public static void SGDMomentum(){
        
    }

    public static void training(int epoch){
        
    }
}
