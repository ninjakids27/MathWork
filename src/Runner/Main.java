package Runner;

import MLComp.*;
import MLComp.ActivationFunctions_Folder.*;
import MLComp.Optimizers_Folder.Optimizers;

public class Main {
    public static void main(String[] args) {
        ColorText.clearScreen();
        Neuron[][] buddy0 = MLOps.buildNeuralNetwork(new int[] { 784, 256, 256, 10 });
        MLOps.heinitilizeWeights(buddy0);
        // overfitting test
        MLOps.training(buddy0,
                "MNIST_CSV/mnist_train.csv",
                0.001,
                15,
                ActivationFunctions::reLU,
                ActivationFunctions::reLUDerivative,
                Optimizers::sgd,
                60_000,
                true);
        MLOps.saveNN(buddy0, "Models//Optimus.ser");

    }
}
