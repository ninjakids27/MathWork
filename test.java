// this is a sandbox where I test stuff. Not actual test cases
public class test{
  public static void main(String[] args){
    int[] a = {5,4,2,2,3};
    Neuron[][] temp = MLOps.buildNeuralNetwork(a);
    double[] weights = {69 ,0.0 ,0.0 ,0.0, 0.0};
    temp[1][0].setWeights(weights);
    MLOps.saveNN(temp, "NN.ser");


    Neuron[][] loadpls = MLOps.loadNN("NN.ser");
    MatrixOps.Print_Vector(loadpls[0][0].getWeights());
  }
}