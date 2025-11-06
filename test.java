// this is a sandbox where I test stuff. Not actual test cases
public class test{
  public static void main(String[] args){
    double[][] A = MatrixOps.generate_Augmented_Matrix(50);

    MatrixOps.RREF(A,true);
  }
}

