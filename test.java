class test{
  public static void main(String[] args){
    double[][] xValues = MatrixOps.generate_Matrix(100, 10);
    double[] yValues = MatrixOps.generate_Random_vector(100);
    double[] coefficients = RegressionFunctions.MLR(xValues, yValues,true);
    // MatrixOps.Print_Vector(coefficients);
  }
}