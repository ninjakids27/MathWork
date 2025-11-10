// this is a sandbox where I test stuff. Not actual test cases
public class test{
  
  public static void main(String[] args){
    double[][] x = MatrixOps.generate_Matrix(10, 11);
    double[] y = MatrixOps.generate_Random_vector(10);


    MatrixOps.Print_Vector(RegressionFunctions.MLR(x,y, true));

  }
}

