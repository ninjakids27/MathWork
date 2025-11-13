// this is a sandbox where I test stuff. Not actual test cases
public class test{
  
  public static void main(String[] args){
    // MLOps.saveNN(MLOps.buildNeuralNetwork(new int[]{784,11,14,16, 10}), "Models/NN784_11_14_16_10.ser");
    // double[] x = {.06,0.135,0.085,0.17,0.285};
    
    // MatrixOps.vectorAdd(x, -0.055);
    // MatrixOps.Print_Vector(x);
    double[] y = {.1,.2,.15,.25,.4};
    RegressionFunctions.stdLinearRegression(x, y);
    double slope = 9.8/RegressionFunctions.linearRegression(x,y)[0];
    System.out.println(slope);
    // System.out.println(Math.pow(slope, -1));

  }
}

