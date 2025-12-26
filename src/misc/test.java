import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

public class test{
  @Test
  public void testDeterminent(){
    double[][] testMatrix = {
      {75,0,0},
      {0,1,0},
      {0,0,1},
    };
    assertEquals(75, MatrixOps.determinant(testMatrix), 0);
  }
  
  @Test
  public void testMultivariableLinearRegression() {
    // Model: y = 3 + 2*x1 - 1*x2
    double[][] X = {
      {0, 0},
      {1, 0},
      {0, 1},
      {2, 3},
      {4, -1}
    };
    double[] y = {
      3,  // 3 + 2*0 -1*0
      5,  // 3 + 2*1 -1*0
      2,  // 3 + 2*0 -1*1
      4,  // 3 + 2*2 -1*3
      12  // 3 + 2*4 -1*(-1)
    };

    // Expect coefficients [intercept, beta1, beta2] = [3, 2, -1]
    // Implement LinearRegression.fit to return this form.
    double[] coeffs = RegressionFunctions.MLR(X, y);
    double[] expected = {3.0, 2.0, -1.0};
    assertArrayEquals(expected, coeffs, 0);
  }


  @Test
  public void testnormalPDF(){
    // https://www.stat.colostate.edu/inmem/gumina/st201/recitation8/downloads/Normal%20Probabilites%20Practice.pdf

    
    // asking for the middle of a normal distrubution should be 0.5 unless floating point bs
    assertEquals(0.5, StatOps.normalPDF(-2, 0, 1),0);
  }

  @Test
  public void testSoftmax(){
        Neuron[][] buddy = MLOps.loadNN("Models//NN784_16_16_10.ser");
        if(buddy == null){
            throw new IllegalArgumentException(ColorText.errorFormat("Could not load model"));
        }
        double[] prob = MLOps.forwardPropagationMNIST(2, "MNIST_CSV/mnist_train.csv",buddy, ActivationFunctions::sigmoid);
        double sum = 0;
        for(double num : prob){
            sum += num;
        }
        assertEquals(1,sum, 0);
  }


  @Test
  public void testLinearRegression(){
    // actual stuff is 11247x+0.113
    double[] x = {1.227,2.12,6.283,4.24,3.347};
    MatrixOps.vectorPow(x, 10e-6);
    double[] y = {0.2,0.4,0.8,0.6,0.5};

    double[] coeffs = RegressionFunctions.linearRegression(x, y);
    assertEquals(11247, coeffs[0], 1);
    assertEquals(0.113, coeffs[1], 0.001);
  }

  @Test
  public void testNeuralNetwork(){
    Neuron[][] buddy = MLOps.loadNN("Models//NN784_128_10.ser");
    if(buddy == null){
        throw new IllegalArgumentException(ColorText.errorFormat("Could not load model"));
    }
    int correct = MLOps.test(buddy, "MNIST_CSV/mnist_test.csv", ActivationFunctions::reLU, ActivationFunctions::reLUDerivative, Optimizers::sgd);
    assertEquals(8900, correct, 100);
  }
}
