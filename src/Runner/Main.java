package Runner;
import MathComp.Matrix;
import MathComp.MatrixOps;
import java.util.HashSet;
import MLComp.Neuron;
import MLComp.MLOps;
public class Main {
    public static void main(String[] args) {
        Neuron[][] network = MLOps.loadNN("Models//Optimus.ser");
        // print the layer number and architecture
        System.out.println("Layer 0: " + network[0].length + " neurons");
        System.out.println("Layer 1: " + network[1].length + " neurons");
        System.out.println("Layer 2: " + network[2].length + " neurons");
        // print the weights of the first neuron in the first layer
    }
}
