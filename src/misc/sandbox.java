public class sandbox {
    public static void main(String[] args){
        Neuron[][] buddy =MLOps.buildNeuralNetwork(new int[]{784,16,16,10});
        MLOps.initilizeWeights(buddy);
        MLOps.saveNN(buddy,"Models/NN784_16_16_10.ser");
    }
}
