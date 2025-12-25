public class sandbox {
    public static void main(String[] args){
        Neuron[][] buddy0 = MLOps.buildNeuralNetwork(new int[] {784,16,16,10});
        MLOps.heinitilizeWeights(buddy0);
        // MatrixOps.printVector(buddy[3][2].getWeights()); 
        MLOps.training(buddy0,
         "MNIST_CSV/mnist_train.csv",
         0.0001,
         50,
         ActivationFunctions::reLU,
         ActivationFunctions::reLUDerivative,
         Optimizers::sgd
        );
    }
}
