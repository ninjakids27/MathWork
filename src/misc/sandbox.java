public class sandbox {
    public static void main(String[] args){
        ColorText.clearScreen();
        Neuron[][] buddy0 = MLOps.buildNeuralNetwork(new int[] {784,128,10});
        MLOps.heinitilizeWeights(buddy0);
        // MatrixOps.printVector(buddy[3][2].getWeights()); 
        MLOps.training(buddy0,
         "MNIST_CSV/mnist_train.csv",
         0.001,
         88.00,
         ActivationFunctions::reLU,
         ActivationFunctions::reLUDerivative,
         Optimizers::adam
        );
        MLOps.saveNN(buddy0, "Models//NN784_16_16_10.ser");
    }
}
