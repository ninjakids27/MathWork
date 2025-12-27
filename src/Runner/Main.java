package Runner;

import MLComp.*;
import MLComp.ActivationFunctions_Folder.*;
import MLComp.Optimizers_Folder.Optimizers;

public class Main {
    public static void main(String[] args) {
        ColorText.clearScreen();
        Neuron[][] buddy0 = MLOps.buildNeuralNetwork(new int[] { 784, 128,64, 10 });
        MLOps.heinitilizeWeights(buddy0);
        // MatrixOps.printVector(buddy[3][2].getWeights());
        MLOps.training(buddy0,
                "MNIST_CSV/mnist_train.csv",
                0.01,
                10,
                ActivationFunctions::reLU,
                ActivationFunctions::reLUDerivative,
                Optimizers::adaGrad,
                60_000,
                true
            );
        MLOps.saveNN(buddy0, "Models//Optimus.ser");

    }
}
