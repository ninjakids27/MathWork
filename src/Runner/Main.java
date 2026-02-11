package Runner;

import MLComp.MLOps;
import MLComp.Neuron;
import MLComp.ActivationFunctions_Folder.ActivationFunctions;
import MLComp.Optimizers_Folder.Optimizers;
public class Main {
    public static void main(String[] args) {
        Neuron[][] testAI = MLOps.loadNN("Models/Optimus.ser");
        int coorect = MLOps.test(
            testAI,
                "MNIST_CSV//mnist_test.csv",
                ActivationFunctions::reLU,
                ActivationFunctions::reLUDerivative
            );

            System.out.println(((double)coorect/100));
    }
}
