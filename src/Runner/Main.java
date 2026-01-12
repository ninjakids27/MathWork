package Runner;

import MLComp.*;
import MLComp.ActivationFunctions_Folder.*;
import MLComp.Optimizers_Folder.Optimizers;
import MathComp.MatrixOps;
public class Main {
    public static void main(String[] args) {
        // ColorText.clearScreen();
        // Neuron[][] buddy0 = MLOps.buildNeuralNetwork(new int[] { 784, 128, 128, 10 });
        // MLOps.heinitilizeWeights(buddy0);
        // // overfitting test
        // MLOps.training(buddy0,
        //         "MNIST_CSV/mnist_train.csv",
        //         0.001,
        //         50,
        //         ActivationFunctions::reLU,
        //         ActivationFunctions::reLUDerivative,
        //         Optimizers::sgd,
        //         60_000
        //     );
        // MLOps.saveNN(buddy0, "Models//Optimus.ser");
        ColorText.clearScreen();
        Neuron[][] buddy0 = MLOps.loadNN("Models//Optimus.ser");
        // double accuracy = (double) MLOps.test(buddy0, "MNIST_CSV/mnist_test.csv", ActivationFunctions::reLU,
        //                 ActivationFunctions::reLUDerivative) / 10000.0;
        //         System.out.println(ColorText.dataFormat("Current Accuracy: ")
        //                 + ColorText.returnFormat(String.format("%.2f", accuracy * 100.0) + "%"));
        MatrixOps.printVector(buddy0[1][0].getWeights());

    }
}
