public class sandbox {
    public static void main(String[] args){
        Neuron[][] buddy = MLOps.loadNN("Models/NN784_16_16_10.ser");
        double[] prob = MLOps.forwardPropagationMNIST(3, "MNIST_CSV/mnist_train.csv",buddy,true);

    }
}
