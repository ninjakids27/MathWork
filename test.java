class test{
  public static void main(String[] args){
    // double[][] x = {{5,8,2,10,3,9,6,1,7,4},{15,10,25,5,20,12,18,30,14,22},{3,5,1,4,6,2,7,8,9,0}};
    // double[] y   = {43,46,51,40,53,44,57,64,63,40};

    // RegressionFunctions.MLR(x, y);

    double[][] a = MatrixOps.generate_Matrix(10, 10);
    MatrixOps.inverseMatrix(a,true);

  }
}