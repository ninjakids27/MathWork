public class sandbox {
    public static void main(String[] args){
        Neuron[][] buddy0 = MLOps.buildNeuralNetwork(new int[] {784,16,16,10});
        MLOps.initilizeWeights(buddy0);
        MLOps.saveNN(buddy0, "Models/NN784_16_16_10.ser");
        Neuron[][] buddy = MLOps.loadNN("Models/NN784_16_16_10.ser");
        // MatrixOps.printVector(buddy[3][2].getWeights()); 
        MLOps.training(buddy0,
         "MNIST_CSV/mnist_train.csv",
         0.0001,
         50,
         ActivationFunctions::leakyReLU,
         ActivationFunctions::leakyReLUDerivative,
         Optimizers::gradientDescent
        );
    }
}
