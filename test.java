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
    assertArrayEquals(expected, coeffs, 1e-9);
  }


}
