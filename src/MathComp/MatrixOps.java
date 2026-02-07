package MathComp;

import java.util.Random;
import Runner.ColorText;
public class MatrixOps {
    private static final double epsilon = 1e-9;
    private static double rounding(double num, int decimal){
        return Math.round(Math.pow(10, decimal)*num)/Math.pow(10, decimal);
    }
    public static double argmax(double[] array){
        double max = array[0];
        int index = 0;
        for(int i = 1; i < array.length; i++){
            if(array[i] > max){
                max = array[i];
                index = i;
            }
        }
        return index;
    }
    /**
     * Generates a random augmented matrix with the specified number of columns.
     * The number of rows will be one less than the number of columns.
     * Each element in the matrix is a random integer between 0 and 9.
     * @param num_columns the number of columns in the augmented matrix
     * @return a 2D array representing the augmented matrix
     */
    public static Matrix<Double> generateAugmentedMatrix(int num_columns) {
        Random random = new Random();
        int num_rows = num_columns - 1;
        Matrix<Double>  Matrix_Name = new Matrix<Double>(num_rows,num_columns);
        int[] coordinates = {0,0};
        for (int k = 0; k < num_rows; k++) {
            coordinates[1] = 0;
            for (int i = 0; i < num_columns; i++) {
                Matrix_Name.set(coordinates, random.nextDouble(10));
                coordinates[1]++;
            }
            coordinates[0]++;
        }
        return Matrix_Name;
    }
    /**
     * Generates a random vector of the specified length.
     * @param length the length of the vector
     * @return a 1D array representing the random vector
     */
    public static double[] generateRandomVector(int length){
        Random random = new Random();
        double[] vector = new double[length];
        for(int i = 0; i < length; i++){
            vector[i] = random.nextInt(10);
        }
        return vector;
    }
    
    /**
     * Accesses a specific row of the given matrix.
     * @param Matrix the input matrix
     * @param row_index the index of the row to access (0-based)
     * @return the specified row as a 1D array
     */
    public static double[] accessRow(double[][] Matrix, int row_index) {
        int num_cols = Matrix[0].length;
        double[] temp = new double[num_cols];
        for (int i = 0; i < num_cols; i++) {
            temp[i] = Matrix[row_index][i];
        }
        return temp;
    }
    /**
     * Accesses a specific column of the given matrix.
     * @param Matrix the input matrix
     * @param column_index the index of the column to access (0-based)
     * @return the specified column as a 1D array
     */
    public static double[] accessColumn(double[][] Matrix, int column_index) {
        int num_rows = Matrix.length;
        double[] temp = new double[num_rows];
        for (int i = 0; i < num_rows; i++) {
            temp[i] = Matrix[i][column_index];
        }
        return temp;
    }
    

    /**
     * Replaces a row in the matrix by adding a multiple of another row to it.
     * @param matrix the input matrix
     * @param replacedRow the index of the row to replace (0-based)
     * @param rowAdded the index of the row to add (0-based)
     * @param factor the factor by which to multiply the added row
     * @return a new matrix with the specified row replaced
     */
    public static Matrix<Double> replacement(Matrix<Double> matrix, int replacedRow, int rowAdded, double factor) {
        Double[] tempRowAdded = matrix.accessRow(rowAdded);
        int numRows = matrix.getRowNum();
        int numCols = matrix.getColumnNum();
        Matrix<Double> tempMatrix = new Matrix<Double>(numRows,numCols);
        // Clone the matrix
        if (factor != 1) {
            constantFactor(tempMatrix, rowAdded, factor);
            tempRowAdded = matrix.accessRow(rowAdded);
        }
        int[] coordinates = {replacedRow,0};
        for (int i = 0; i < tempRowAdded.length; i++) {
            tempMatrix.set(coordinates, matrix.get(coordinates)+tempRowAdded[i]);
            // tempMatrix[replacedRow][i] = matrix[replacedRow][i] + tempRowAdded[i];
            coordinates[1]++;
        }
        return tempMatrix;
    }

    /**
     * Multiplies a scalar with a vector.
     * @param a the scalar value
     * @param vector the input vector
     * @return the resulting vector after multiplication
     */
    public static double[] scalarVector(double a,double[] vector){
		double[] tempVector = new double[vector.length];
		for(int i = 0; i < vector.length;i++){
			tempVector[i] = a*vector[i];
		}
		return tempVector;
	}

    /**
     * Multiplies a scalar with a matrix.
     * @param a the scalar value
     * @param matrix the input matrix
     * @return the resulting matrix after multiplication
     */
	public static double[][] scalarMatrix(double a, double[][] matrix){
		double[][] tempMatrix = new double[matrix.length][matrix[0].length];
		for(int i = 0; i < matrix.length;i++){
			for(int k = 0; k < matrix[0].length;k++){
				tempMatrix[i][k] = a*matrix[i][k];
			}
		}
		return tempMatrix;
	}

    /**
     * Multiplies a specific row of the matrix by a constant factor and returns a new matrix.
     * @param Matrix the input matrix
     * @param row_index the index of the row to modify (0-based)
     * @param factor the factor by which to multiply the row
     * @return a new matrix with the specified row modified
     */
    public static void constantFactor(Matrix<Double> Matrix, int rowIndex, double factor) {
        Double[] temp = Matrix.accessRow(rowIndex);
        for(int i = 0 ; i < temp.length; i++){
            temp[i]*=factor;
        }
        Matrix.setRow(rowIndex, temp);
    }

    /**
     * Prints the given vector to the console.
     * @param Vector the input vector
     */
    public static void printVector(double[] Vector) {
        for (int i = 0; i < Vector.length; i++) {
            System.out.print(Vector[i] + " ");
        }
        System.out.println("");
    }
    public static void printVector(int[] Vector) {
        for (int i = 0; i < Vector.length; i++) {
            System.out.print(Vector[i] + " ");
        }
        System.out.println("");
    }

    /**
     * Determines the type of solution for a system of linear equations represented by the given RREF matrix.
     * @param RREFMatrix the RREF matrix
     * @return a string indicating the type of solution ("No Solution", "Unique Solution", or "Infinitely Many Solutions (Free Variables)")
     */
    public static String getSolutionType(Matrix<Double> RREFMatrix) {
    int numRows = RREFMatrix.getRowNum();
    int numCols = RREFMatrix.getColumnNum();
    int numVariables = numCols - 1; // Assuming it's an augmented matrix

    // Check for no solution
    int[] coordinates = {0,0};
    for (int i = 0; i < numRows; i++) {
        coordinates[0] = i;
        boolean allZerosExceptLast = true;
        for (int j = 0; j < numVariables; j++) {
            coordinates[1] = j;
            if (RREFMatrix.get(coordinates) != 0.0) {
                allZerosExceptLast = false;
                break;
            }
        }
        if (allZerosExceptLast && RREFMatrix.get(new int[] {i,numVariables}) != 0.0) {
            return "No Solution";
        }
    }

    // Check for free variables
    coordinates = new int[] {0,0};
    int pivotCount = 0;
    for (int i = 0; i < numRows; i++) {
        coordinates[0] = i;
        for (int j = 0; j < numVariables; j++) {
            coordinates[1] = j;
            if (RREFMatrix.get(coordinates) == 1.0) {
                // Check if it's a leading 1
                boolean isLeadingOne = true;
                for (int k = 0; k < j; k++) {
                    coordinates[1] = k;
                    if (RREFMatrix.get(coordinates) != 0.0) {
                        isLeadingOne = false;
                        break;
                    }
                }
                if (isLeadingOne) {
                    pivotCount++;
                    break;
                }
            }
        }
    }

    if (pivotCount < numVariables) {
        return "Infinitely Many Solutions (Free Variables)";
    }

    return "Unique Solution";
    }
    /**
     * Converts the given matrix to its Reduced Row Echelon Form (RREF).
     * @param Matrix the input matrix
     * @param debug if true, prints intermediate steps for debugging
     * @return the RREF of the input matrix
     */
    public static Matrix<Double> RREF(Matrix<Double> Matrix, boolean debug) {
    int numRows = Matrix.getRowNum();
    int numCols = Matrix.getColumnNum();
    Matrix<Double> tempMatrix= new Matrix<Double>(numRows,numCols);
    int[] coordinates = {0,0};
    for (int i = 0; i < numRows; i++) {
        coordinates[0] = i;
        for (int k = 0; k < numCols; k++) {
            coordinates[1] = k;
            tempMatrix.set(coordinates, Matrix.get(coordinates));
            
        }
    }


    int pivotRow = 0;
    coordinates = new int[] {0,0};
    for (int j = 0; j < numCols && pivotRow < numRows; j++) {
        coordinates[1] = j;
        int maxRow = pivotRow;
        for (int i = pivotRow + 1; i < numRows; i++) {
            coordinates[0] = i;
            if (Math.abs(tempMatrix.get(coordinates)) > Math.abs(tempMatrix.get(new int[] {maxRow, coordinates[1]}))) {
                maxRow = i;
            }
        }

        if (Math.abs(tempMatrix.get(new int[] {maxRow, coordinates[1]})) > epsilon) {
            tempMatrix.switchRows(pivotRow, maxRow);
            if(debug){
                therefore();
                System.out.println(tempMatrix);
            }

            double pivotValue = tempMatrix.get(new int[] {pivotRow, j});
            constantFactor(tempMatrix, pivotRow, 1/pivotValue);
            roundMatrix(tempMatrix, epsilon);
            if(debug){
                therefore();
                System.out.println(tempMatrix);
            }

            for(int i = pivotRow + 1; i < numRows; i++){
                double factor = tempMatrix.get(new int[] {i, j});
                replacement(tempMatrix, i, pivotRow, -factor);
                roundMatrix(tempMatrix, epsilon);
                if(debug){
                    therefore();
                    System.out.println(tempMatrix);
                }
            }
            pivotRow++;
        }
    }
    coordinates = new int[] {0,0};
    for(int i = numRows - 1; i >=0; i--){
        int pivotCol = -1;
        coordinates[0] = i;
        for(int j = 0; j < numCols; j++){
            coordinates[1] = j;
            if(Math.abs(tempMatrix.get(coordinates) - 1.0) < epsilon){ 
                pivotCol = j;
                break;
            }
        }
        if(pivotCol != -1){
            for(int k = 0; k < i; k++){
                coordinates[0] = k;
                coordinates[1] = pivotCol;
                double factor = tempMatrix.get(coordinates);
                replacement(tempMatrix, k, i, -factor);
                roundMatrix(tempMatrix, epsilon);
                if(debug){
                    therefore();
                    System.out.println(tempMatrix);
                }
            }
        }
    }
    
    roundMatrix(tempMatrix, epsilon);
    if(debug){
        therefore();
    }
    return tempMatrix;
}
    public static Matrix<Double> RREF(Matrix<Double> Matrix) {
    int numRows = Matrix.getRowNum();
    int numCols = Matrix.getColumnNum();
    Matrix<Double> tempMatrix = new Matrix<Double>(numRows,numCols);
    int[] coordinates = {0,0};
    for (int i = 0; i < numRows; i++) {
        coordinates[0] = i;
        for (int k = 0; k < numCols; k++) {
            coordinates[1] = k;
            tempMatrix.set(coordinates, Matrix.get(coordinates));
        }
    }

    
    int pivotRow = 0;
    coordinates = new int[] {0,0};
    for (int j = 0; j < numCols && pivotRow < numRows; j++) {
        int maxRow = pivotRow;
        coordinates[1] = j;
        for (int i = pivotRow + 1; i < numRows; i++) {
            coordinates[0] = i;
            if (Math.abs(tempMatrix.get(coordinates)) > Math.abs(tempMatrix.get(new int[] {maxRow, coordinates[1]}))) {
                maxRow = i;
            }
        }

        if (Math.abs(tempMatrix.get(new int[] {maxRow, coordinates[1]})) > epsilon) { 
            tempMatrix.switchRows(pivotRow, maxRow);

            double pivotValue = tempMatrix.get(new int[] {pivotRow, j});
            constantFactor(tempMatrix, pivotRow, 1/pivotValue);
            roundMatrix(tempMatrix, epsilon);

            for(int i = pivotRow + 1; i < numRows; i++){
                double factor = tempMatrix.get(new int[] {i, j});
                replacement(tempMatrix, i, pivotRow, -factor);
                roundMatrix(tempMatrix, epsilon);
            }
            pivotRow++;
        }
    }
    
    for(int i = numRows - 1; i >=0; i--){
        int pivotCol = -1;
        for(int j = 0; j < numCols; j++){
            if(Math.abs(tempMatrix.get(new int[] {i, j}) - 1.0) < epsilon){ 
                pivotCol = j;
                break;
            }
        }
        if(pivotCol != -1){
            for(int k = 0; k < i; k++){
                double factor = tempMatrix.get(new int[] {k, pivotCol});
                replacement(tempMatrix, k, i, -factor);
                roundMatrix(tempMatrix, epsilon);
            }
        }
    }
    
    roundMatrix(tempMatrix, epsilon);
    return tempMatrix;
}
    /**
     * Rounds elements of the matrix that are smaller than epsilon to zero.
     * @param matrix the input matrix
     * @param epsilon the threshold below which values are rounded to zero
     * @return a new matrix with small values rounded to zero
     */
    private static void roundMatrix(Matrix<Double> matrix, double epsilon) {
    int numRows = matrix.getRowNum();
    int numCols = matrix.getColumnNum();
    int[] coordinates = {0,0};
    for (int i = 0; i < numRows; i++) {
        coordinates[0] = i;
        for (int j = 0; j < numCols; j++) {
            coordinates[1] = j;
            if (Math.abs(matrix.get(coordinates)) < epsilon) {
                matrix.set(coordinates, 0.0);
            }
        }
    }
}   
    /**
     * Prints an arrow symbol ("->") to the console.
     */
    public static void therefore() {
        System.out.println("->");
    }
    /**
     * Adds two vectors element-wise.
     * @param vector1 the first vector
     * @param vector2 the second vector
     * @return the element-wise sum of the two vectors
     */
    public static double[] vectorSum(double[] vector1, double[] vector2){
        
        if(vector1.length != vector2.length){
            throw new IllegalArgumentException(ColorText.errorFormat("Mismatched lengths"));
        }

        double[] tempVector = new double[vector1.length];
        for(int i = 0; i < vector1.length; i++){
            tempVector[i] = vector1[i]+vector2[i];
        }
        return tempVector;
    }
    /**
     * Multiplies a matrix by a vector.
     * @param matrix the input matrix
     * @param vector the input vector
     * @return the resulting vector after multiplication
     */
    public static double[] matrixVectorMult(Matrix<Double> matrix, double[] vector){
        if(matrix.getColumnNum() != vector.length){
            throw new IllegalArgumentException(ColorText.errorFormat("Incompatible matrix and vector dimensions"));
        }

        double[] result = new double[matrix.getRowNum()];
        for (int i = 0; i < matrix.getRowNum(); i++) {
            result[i] = dotProduct(matrix.accessRow(i), vector);
        }
        return result;
    }
    /**
     * Calculates the dot product of two vectors.
     * @param vector1 the first vector
     * @param vector2 the second vector
     * @param debug if true, prints intermediate products
     * @return the dot product of the two vectors
     */
    public static double dotProduct(double[] vector1,double[] vector2, boolean debug){
        double temp = 0;
        if(vector1.length == vector2.length){
            for(int i = 0; i < vector1.length; i++){
                temp += vector1[i]*vector2[i];
                if(debug){
                    System.out.println(vector1[i]*vector2[i]);
                }
            }
            return temp;
        }else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input vectors are wrong lengths"));
        }
    }
    /**
     * Converts nanoseconds to seconds.
     * @param nanoseconds the time in nanoseconds
     * @return the time in seconds
     */
    
    public static double nanoToSeconds(long nanoseconds){
        return nanoseconds/1e9;
    }
    
    public static double dotProduct(double[] vector1,double[] vector2){
        double temp = 0;
        if(vector1.length == vector2.length){
            for(int i = 0; i < vector1.length; i++){
                temp += vector1[i]*vector2[i];
            }
            return temp;
        }else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input vectors are wrong lengths"));
        }
    }

    public static double dotProduct(Double[] vector1,double[] vector2){
        double temp = 0;
        if(vector1.length == vector2.length){
            for(int i = 0; i < vector1.length; i++){
                temp += vector1[i]*vector2[i];
            }
            return temp;
        }else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input vectors are wrong lengths"));
        }
    }

    public static double dotProduct(Double[] vector1,Double[] vector2){
        double temp = 0;
        if(vector1.length == vector2.length){
            for(int i = 0; i < vector1.length; i++){
                temp += vector1[i]*vector2[i];
            }
            return temp;
        }else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input vectors are wrong lengths"));
        }
    }
    
    public static int dotProduct(int[] vector1,int[] vector2){
        int temp = 0;
        if(vector1.length == vector2.length){
            for(int i = 0; i < vector1.length; i++){
                temp += vector1[i]*vector2[i];
            }
            return temp;
        }else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input vectors are wrong lengths"));
        }
    }

    /**
     * Generates a random matrix with the specified number of rows and columns.
     * Each element in the matrix is a random integer between 0 and 9.
     * @param Rows the number of rows in the matrix
     * @param Columns the number of columns in the matrix
     * @return a 2D array representing the random matrix
     */
    public static double[][] generateMatrix(int Rows, int Columns){
        double[][] Matrix = new double[Rows][Columns];
        Random random = new Random();
        for(int i = 0; i < Rows; i++){
            for(int k = 0; k < Columns; k++){
                Matrix[i][k] = random.nextInt(10);
            }
        }
        return Matrix;

    }
    
    /**
     * Multiplies two matrices.
     * @param MatrixA the first matrix
     * @param MatrixB the second matrix
     * @param debug if true, prints intermediate dot products
     * @return the resulting matrix after multiplication
     */
    public static double[][] matrixMult(double[][] MatrixA, double[][] MatrixB, boolean debug){
        if(MatrixA[0].length == MatrixB.length){
            double[][] ResultMatrix = new double[MatrixA.length][MatrixB[0].length];
            for(int i = 0; i < MatrixA.length; i++){
                for(int j = 0; j < MatrixB[0].length; j++){
                    double[] Row = accessRow(MatrixA, i);
                    double[] Column = accessColumn(MatrixB, j);
                    ResultMatrix[i][j] = dotProduct(Row, Column, debug);
                }
            }
            return ResultMatrix;
        } else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input Matrices are invalid"));

        }
    }

    /**
     * Multiplies two matrices. O(n^3) time complexity.
     * @param Matrix_A
     * @param Matrix_B
     * @return
     */
    public static Matrix<Double> matrixMult(Matrix<Double> MatrixA, Matrix<Double> MatrixB){
        if(MatrixA.getColumnNum() == MatrixB.getRowNum()){
            int[] coordinates = {0,0};
            Matrix<Double> ResultMatrix = new Matrix<Double>(MatrixA.getRowNum(), MatrixB.getColumnNum());
            for(int i = 0; i < MatrixA.getRowNum(); i++){
                coordinates[0] = i;
                for(int j = 0; j < MatrixB.getColumnNum(); j++){
                    coordinates[1] = j;
                    Double[] Row = MatrixA.accessRow(i);
                    Double[] Column = MatrixB.accessColumn(j);
                    ResultMatrix.set(coordinates, dotProduct(Row, Column));
                }
            }
            return (ResultMatrix);
        } else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input Matrices are invalid"));
        }
    }

    /**
     * Calculates the cross product of two 3D vectors.
     * @param vector1 the first vector (length 3)
     * @param vector2 the second vector (length 3)
     * @return the resulting vector after cross product
     */
    public static double[] crossProduct(double[] vector1, double[] vector2){
        if(vector1.length != 3 || vector2.length != 3)
            throw new IllegalArgumentException("");

        double[] result_vector = {
            (vector1[2]*vector2[1]-vector1[1]*vector2[2]),
            (vector1[0]*vector2[2]-vector1[2]*vector2[0]),
            (vector1[1]*vector2[0]-vector1[0]*vector2[1])
        };
        return result_vector;
    }

    /**
     * Transposes the given matrix.
     * @param matrix the input matrix
     * @return the transposed matrix
     */
    public static Matrix<Double> matrixTransposition(Matrix<Double> matrix){
        Matrix<Double> tempMatrix = new Matrix<Double>(matrix.getColumnNum(), matrix.getRowNum());
        Double[] tempRow = new Double[tempMatrix.getColumnNum()];
        int[] coordinates = {0,0};
        for(int i = 0; i<matrix.getRowNum();i++){
            coordinates[0] = i;
            tempRow = matrix.accessRow(i);
            for(int j = 0; j < tempMatrix.getRowNum(); j++){
                coordinates[1] = j;
                tempMatrix.set(coordinates, tempRow[j]);
            }
        }
        return tempMatrix;
    }

    /**
     * Generates a submatrix by removing the first row and the specified column from the given matrix.
     * @param matrix the input matrix
     * @param column the column to remove
     * @return the resulting submatrix
     */
    public static Matrix<Double> subMatrix(Matrix<Double> matrix, int column) {
        Matrix<Double> tempMatrix = new Matrix<Double>(matrix.getRowNum() - 1, matrix.getColumnNum() - 1);
        int[] coordinates = {0, 0};
        for (int i = 1; i < matrix.getRowNum(); i++) {
            coordinates[0] = i - 1;
            int subCol = 0;
            for (int j = 0; j < matrix.getColumnNum(); j++) {
                coordinates[1] = subCol;
                if (j == column) continue;
                tempMatrix.set(coordinates, matrix.get(new int[] {i, j}));
                subCol++;
            }
        }
        return tempMatrix;
    }

    /**
     * Calculates the determinant of a square matrix using recursion.
     * @param matrix the input square matrix
     * @return the determinant of the matrix
     */
    public static double determinant(Matrix<Double> matrix){
		double sum = 0;
		if(matrix.getRowNum() != matrix.getColumnNum()){
			throw new IllegalArgumentException(ColorText.errorFormat("determinent input is not equal"));
		}
		int n = matrix.getRowNum();
		
		if(n==1){
			return matrix.get(new int[] {0, 0});
		}
		if(n == 2){
			return (matrix.get(new int[] {0, 0}) * matrix.get(new int[] {1, 1})-
			matrix.get(new int[] {0, 1}) * matrix.get(new int[] {1, 0}));
		}
        double sign;
        for(int i = 0; i < n; i++){
            if(i%2==0){
                sign = 1.0;
            }else{
                sign = -1.0;
            }
            sum += sign* matrix.get(new int[] {0, i}) * determinant(subMatrix(matrix, i));
        }
        return sum;
	}

    /**
     * Calculates the inverse of a square matrix using Gaussian elimination.
     * @param Matrix the input square matrix
     * @param stopwatch if true, measures and prints the time taken for the operation
     * @return the inverse of the matrix
     */
    public static Matrix<Double> inverseMatrix(Matrix<Double> Matrix){
        if(Matrix.getRowNum() != Matrix.getColumnNum()){
            throw new IllegalArgumentException(ColorText.errorFormat("Not a square matrix"));
        }
        if(determinant(Matrix) == 0){
            throw new IllegalArgumentException(ColorText.errorFormat("Determinent is zero"));
        }
        int n = Matrix.getRowNum();

        if(n == 2){
            double det = determinant(Matrix);
            Double[] data = new Double[] {
                Matrix.get(new int[] {1, 1})/det, -Matrix.get(new int[] {0, 1})/det,
                -Matrix.get(new int[] {1, 0})/det, Matrix.get(new int[] {0, 0})/det
            };
            Matrix<Double> inverse = new Matrix<Double>(2,2, data);
            return inverse;
        }
        int numRows = Matrix.getRowNum();
        int numCols = Matrix.getColumnNum();
        Matrix<Double> tempMatrix = new Matrix<Double>(numRows,numCols);
        Matrix<Double> inverseMatrix = new Matrix<Double>(numRows,numCols);
        for(int i = 0; i < numRows; i++) {
            for (int k = 0; k < numCols; k++) {
                tempMatrix.set(new int[] {i, k}, Matrix.get(new int[] {i, k}));
            }
        }

        for(int i = 0; i < numRows;i++){
            inverseMatrix.set(new int[] {i, i}, 1.0);
        }


        
        int pivotRow = 0;

        for (int j = 0; j < numCols && pivotRow < numRows; j++) {
            int maxRow = pivotRow;
            for (int i = pivotRow + 1; i < numRows; i++) {
                if (Math.abs(tempMatrix.get(new int[] {i, j})) > Math.abs(tempMatrix.get(new int[] {maxRow, j}))) {
                    maxRow = i;
                }
            }

            if (Math.abs(tempMatrix.get(new int[] {maxRow, j})) > epsilon) { 
                tempMatrix.switchRows(pivotRow, maxRow);
                inverseMatrix.switchRows(pivotRow, maxRow);

                double pivotValue = tempMatrix.get(new int[] {pivotRow, j});
                constantFactor(tempMatrix, pivotRow, 1/pivotValue);
                roundMatrix(tempMatrix, epsilon);
                constantFactor(inverseMatrix, pivotRow, 1/pivotValue);
                roundMatrix(inverseMatrix, epsilon);

                for(int i = pivotRow + 1; i < numRows; i++){
                    double factor = tempMatrix.get(new int[] {i, j});
                    replacement(tempMatrix, i, pivotRow, -factor);
                    roundMatrix(tempMatrix, epsilon);
                    replacement(inverseMatrix, i, pivotRow, -factor);
                    roundMatrix(inverseMatrix, epsilon);
                }
                pivotRow++;
            }
        }
        
        for(int i = numRows - 1; i >=0; i--){
            int pivotCol = -1;
            for(int j = 0; j < numCols; j++){
                if(Math.abs(tempMatrix.get(new int[] {i, j}) - 1.0) < epsilon){ 
                    pivotCol = j;
                    break;
                }
            }
            if(pivotCol != -1){
                for(int k = 0; k < i; k++){
                    double factor = tempMatrix.get(new int[] {k, pivotCol});
                    replacement(tempMatrix, k, i, -factor);
                    roundMatrix(tempMatrix, epsilon);
                    replacement(inverseMatrix, k, i, -factor);
                    roundMatrix(inverseMatrix, epsilon);
                }
            }
        }
        
        roundMatrix(inverseMatrix, epsilon);
        return inverseMatrix;
    }


    public static Matrix<Double> inverseMatrix(Matrix<Double> Matrix,boolean stopwatch){
        long startTime = System.nanoTime();
        
        Matrix<Double> inverseMatrix = inverseMatrix(Matrix);
        if(stopwatch){
            long estimatedTime = System.nanoTime() - startTime;
            System.out.println("elapsed time in seconds: "+nanoToSeconds(estimatedTime));
        }
        return inverseMatrix;
    }
    // x vector returns y vector where y1 = x1^factor
    /**
     * Raises each element in the input vector to the specified power.
     * @param vector the data vector
     * @param factor the power to raise each element to
     * @return the transformed vector
     */
    public static double[] vectorPow(double[] vector, double factor){
        for(int i = 0; i<vector.length; i++){
            vector[i] = Math.pow(vector[i], factor);
        }
        return vector;
    }
    
    /**
     * Adds a specified factor to each element in the input vector.
     * @param vector the data vector
     * @param factor the factor to add to each element
     * @return the transformed vector
     */
    public static double[] vectorAdd(double[] vector, double factor){
        for(int i = 0; i<vector.length; i++){
            vector[i] += factor;
        }
        return vector;
    }

    /**
     * ln(vector) ;)
     * @param vector double[] type
     * @return  double of log
     */
    public static double[] vectorLog(double[] vector){
        double temp[] = new double[vector.length];
        for(int i = 0; i<vector.length;i++){
            temp[i] = Math.log(vector[i]);
        }
        return temp;
    }

    /**
     * Multiplies each element in the input vector by the specified factor.
     * @param vector the data vector
     * @param factor the factor to multiply each element by
     * @return the transformed vector
     */
    public static double[] vectorMult(double[] vector, double factor){
        for(int i = 0; i<vector.length; i++){
            vector[i] *= factor;
        }
        return vector;
    }
    
    public static double[] vectorMult( double factor,double[] vector){
        for(int i = 0; i<vector.length; i++){
            vector[i] *= factor;
        }
        return vector;
    }
    /**
     * Calculates the element-wise product of two vectors.
     * @param vector1 the first data vector
     * @param vector2 the second data vector
     * @return the element-wise product as a new vector
     */
    public static double[] vectorProduct(double[] vector1, double[] vector2){
        double[] tempVector = new double[vector1.length];
        if(vector1.length != vector2.length){
            throw new IllegalArgumentException(ColorText.errorFormat("Fucked up lengths"));
        }
        for(int i = 0; i<vector1.length; i++){
            tempVector[i] = vector1[i]*vector2[i];
        }
        return tempVector;
    }


    /**
     * Rounds each element in the input vector to the specified number of decimal places.
     * @param num the input vector
     * @param decimal the number of decimal places to round to
     * @return the rounded vector
     */
    public static double[] roundingVector(double[] num, int decimal){
        double[] roundedVector = new double[num.length];
            for(int i = 0; i < num.length; i++){
                roundedVector[i] = rounding(num[i],decimal);
            }
            return roundedVector;
    }

    public static Matrix<Double> hadamardProduct(Matrix<Double> A, Matrix<Double> B){
        if(A.getRowNum() != B.getRowNum() || A.getColumnNum() != B.getColumnNum())
            throw new IllegalArgumentException(ColorText.errorFormat("Matrices are not equal length"));
        
        for(int i = 0; i < A.getRowNum(); i++){
            for(int j = 0; j < A.getColumnNum(); j++){
                A.set(i, j, A.get(i, j) * B.get(i, j));
            }
        }
        return A;
    }
}