class test{
  // not an actual test file just a sandbox to test code snippets
  public static void main(String[] args){
      double[] x = {0.14,0.26,0.64,0.765,0.89};
      double[] y = {0.31,0.4,0.69,0.74,0.8};
      RegressionFunctions.stdPowerRegression(x,y);
      MatrixOps.Print_Vector(RegressionFunctions.powerRegression(x,y));
  }
}