public class sandbox {
    public static void main(String[] args){
        Neuron[][] buddy = MLOps.loadNN("C:\\githubProject\\MathWork\\Models\\NN784_16_16_10.ser");
        MLOps.initilizeWeights(buddy);
        MatrixOps.Print_Vector(buddy[2][0].getWeights());
    }
}
