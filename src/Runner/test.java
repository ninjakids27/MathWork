package Runner;

import MathComp.*;
import MLComp.*;
import MLComp.ActivationFunctions_Folder.*;
import MLComp.Optimizers_Folder.*;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

public class test {
  @Test
  public void testDeterminent() {
    double[][] testMatrix = {
        { 75, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
    };
    assertEquals(75, MatrixOps.determinant(testMatrix), 0);
  }

  @Test
  public void testMultivariableLinearRegression() {
    // Model: y = 3 + 2*x1 - 1*x2
    double[][] X = {
        { 0, 0 },
        { 1, 0 },
        { 0, 1 },
        { 2, 3 },
        { 4, -1 }
    };
    double[] y = {
        3, // 3 + 2*0 -1*0
        5, // 3 + 2*1 -1*0
        2, // 3 + 2*0 -1*1
        4, // 3 + 2*2 -1*3
        12 // 3 + 2*4 -1*(-1)
    };

    // Expect coefficients [intercept, beta1, beta2] = [3, 2, -1]
    // Implement LinearRegression.fit to return this form.
    double[] coeffs = RegressionFunctions.MLR(X, y);
    double[] expected = { 3.0, 2.0, -1.0 };
    assertArrayEquals(expected, coeffs, 0);
  }

  @Test
  public void testnormalPDF() {
    // https://www.stat.colostate.edu/inmem/gumina/st201/recitation8/downloads/Normal%20Probabilites%20Practice.pdf

    // asking for the middle of a normal distrubution should be 0.5 unless floating
    // point bs
    assertEquals(0.5, StatOps.normalPDF(-2, 0, 1), 0);
  }

  @Test
  public void testSoftmax() {
    Neuron[][] buddy = MLOps.loadNN("Models//NN784_16_16_10.ser");
    if (buddy == null) {
      throw new IllegalArgumentException(ColorText.errorFormat("Could not load model"));
    }
    double[] prob = MLOps.forwardPropagationMNIST(2, "MNIST_CSV/mnist_train.csv", buddy, ActivationFunctions::sigmoid);
    double sum = 0;
    for (double num : prob) {
      sum += num;
    }
    assertEquals(1, sum, 0);
  }

  @Test
  public void testLinearRegression() {
    // actual stuff is 11247x+0.113
    double[] x = { 1.227, 2.12, 6.283, 4.24, 3.347 };
    MatrixOps.vectorPow(x, 10e-6);
    double[] y = { 0.2, 0.4, 0.8, 0.6, 0.5 };

    double[] coeffs = RegressionFunctions.linearRegression(x, y);
    assertEquals(11247, coeffs[0], 1);
    assertEquals(0.113, coeffs[1], 0.001);
  }

  @Test
  public void testNeuralNetwork() {
    Neuron[][] buddy = MLOps.loadNN("Models//Optimus.ser");
    if (buddy == null) {
      throw new IllegalArgumentException(ColorText.errorFormat("Could not load model"));
    }
    int correct = MLOps.test(buddy, "MNIST_CSV/mnist_test.csv", ActivationFunctions::reLU,
        ActivationFunctions::reLUDerivative);
    assertEquals(8800, correct, 100);
  }


 @Test
 public void testTensorGetCommand(){
   Tensor<Integer> a = new Tensor<Integer>(new Integer[] {1,2,3,4,5,6,7,8}, new int[] {2,2,2});
   int[] coordinates = {0,0,0};
   int result  = a.get(coordinates);
   assertEquals(result,  1, 0);
   coordinates = new int[] {2,2,2};
   result  = a.get(coordinates);
   assertEquals(result,  8, 0);
 }


 @Test
    public void testTensorIndexing() {
        // 1. Define dimensions: 2 sheets, 3 rows, 4 columns (2x3x4)
        int[] dims = {2,2,2,2,5,5,5};
        
        // 2. Create the flat data array (2 * 3 * 4 = 24 elements)
        Integer[] flatData = new Integer[2000];
        for (int i = 0; i < 2000; i++) {
            flatData[i] = i; // Fill with an arbitrary large array
        }

        // 3. Initialize your Tensor class
        Tensor<Integer> tensor = new Tensor<>(flatData, dims);

        // test all get cases
        for(int i = 0; i < 2000; i++){
          assertEquals((Integer)i, tensor.get(new int[] {
            (i / (1000)) % 2,
            (i / (500)) % 2,
            (i / (250)) % 2,
            (i / (125)) % 2,
            (i / (25)) % 5,
            (i / (5)) % 5,
            i % 5
          }));
          tensor.set(new int[] {
            (i / (1000)) % 2,
            (i / (500)) % 2,
            (i / (250)) % 2,
            (i / (125)) % 2,
            (i / (25)) % 5,
            (i / (5)) % 5,
            i % 5,
          },1);
        }
        
        // test all set cases
        for(int i = 0; i < 2000; i++){
          assertEquals((Integer)1, tensor.get(new int[] {
            (i / (1000)) % 2,
            (i / (500)) % 2,
            (i / (250)) % 2,
            (i / (125)) % 2,
            (i / (25)) % 5,
            (i / (5)) % 5,
            i % 5
          }));
        }
    }

  
}
