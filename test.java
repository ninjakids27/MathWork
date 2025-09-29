class test{
  public static void main(String[] args){
    // Independent variables (features)
        double[][] independentVariables = {
            {2.5, 5.0},
            {3.1, 5.8},
            {4.2, 7.1},
            {5.5, 9.3},
            {6.8, 11.2},
            {7.9, 13.0},
            {8.5, 14.5},
            {9.0, 15.1},
            {10.1, 16.8},
            {11.5, 18.5}
        };

        // Dependent variable (target)
        double[] dependentVariable = {
            12.6,
            15.5,
            19.8,
            24.1,
            29.5,
            34.0,
            38.5,
            40.2,
            45.1,
            49.8
        };
    // System.out.println(MatrixOps.determinant(independentVariables));
    MatrixOps.Print_Vector(RegressionFunctions.MLR(independentVariables, dependentVariable));

    // double[][] a = MatrixOps.generate_Matrix(10, 10);
    // MatrixOps.inverseMatrix(a,true);






        
    // double[] a = {-70_000,250_000};
    // double[] b = {0.6,0.4};
    // StatOps.probability_distribution(a,b,2);
  }
}