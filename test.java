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
  public void testMean(){
    double[] a = {5,3,6,7,8};
    assertEquals(5.8,StatOps.mean(a),0);
  }
}
