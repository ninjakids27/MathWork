import java.util.Random;

public class MatrixOps {
    private static final double epsilon = 1e-9; 
    public static double[][] generate_Augmented_Matrix(int num_columns) {
        Random random = new Random();
        int num_rows = num_columns - 1;
        double[][] Matrix_Name = new double[num_rows][num_columns];
        for (int k = 0; k < num_rows; k++) {
            for (int i = 0; i < num_columns; i++) {
                Matrix_Name[k][i] = random.nextInt(10);
            }
        }
        return Matrix_Name;
    }
    public static double[] generate_Random_vector(int length){
        Random random = new Random();
        double[] vector = new double[length];
        for(int i = 0; i < length; i++){
            vector[i] = random.nextInt(10);
        }
        return vector;
    }
    public static void Print_Matrix(double[][] Matrix_Name) {
        int num_rows = Matrix_Name.length;
        int num_columns = Matrix_Name[0].length;
        for (int i = 0; i < num_rows; i++) {
            for (int k = 0; k < num_columns; k++) {
                System.out.print(Matrix_Name[i][k] + " ");
            }
            System.out.println("");
        }
    }
    public static double[] access_row(double[][] Matrix, int row_index) {
        int num_cols = Matrix[0].length;
        double[] temp = new double[num_cols];
        for (int i = 0; i < num_cols; i++) {
            temp[i] = Matrix[row_index][i];
        }
        return temp;
    }
    public static double[] access_column(double[][] Matrix, int column_index) {
        int num_rows = Matrix.length;
        double[] temp = new double[num_rows];
        for (int i = 0; i < num_rows; i++) {
            temp[i] = Matrix[i][column_index];
        }
        return temp;
    }
    public static double[][] switch_rows(double[][] Matrix, int row1_index, int row2_index) {
        // Clone the matrix
        int numRows = Matrix.length;
        int numCols = Matrix[0].length;
        double[][] temp_Matrix = new double[numRows][numCols];
        for (int i = 0; i < numRows; i++) {
            for (int k = 0; k < numCols; k++) {
                temp_Matrix[i][k] = Matrix[i][k];
            }
        }
        double[] row1 = access_row(temp_Matrix, row1_index);
        double[] row2 = access_row(temp_Matrix, row2_index);
        for (int i = 0; i < row1.length; i++) {
            temp_Matrix[row1_index][i] = row2[i];
            temp_Matrix[row2_index][i] = row1[i];
        }
        return temp_Matrix;
    }
    public static double[][] replacement(double[][] matrix, int replacedRow, int rowAdded, double factor) {
        double[] tempRowAdded = access_row(matrix, rowAdded);
        int numRows = matrix.length;
        int numCols = matrix[0].length;
        double[][] tempMatrix = new double[numRows][numCols];
        // Clone the matrix
        for (int i = 0; i < numRows; i++) {
            for (int k = 0; k < numCols; k++) {
                tempMatrix[i][k] = matrix[i][k];
            }
        }
        if (factor != 1) {
            tempRowAdded = access_row(Constant_factor(matrix, rowAdded, factor), rowAdded);
        }
        for (int i = 0; i < tempRowAdded.length; i++) {
            tempMatrix[replacedRow][i] = matrix[replacedRow][i] + tempRowAdded[i];
        }
        return tempMatrix;
    }
    	public static double[] scalarVector(double a,double[] vector){
		double[] tempVector = new double[vector.length];
		for(int i = 0; i < vector.length;i++){
			tempVector[i] = a*vector[i];
		}
		return tempVector;
	}
	public static double[][] scalarMatrix(double a, double[][] matrix){
		double[][] tempMatrix = new double[matrix.length][matrix[0].length];
		for(int i = 0; i < matrix.length;i++){
			for(int k = 0; k < matrix[0].length;k++){
				tempMatrix[i][k] = a*matrix[i][k];
			}
		}
		return tempMatrix;
	}

    public static double[][] Constant_factor(double[][] Matrix, int row_index, double factor) {
        double[] temp = access_row(Matrix, row_index);
        int numRows = Matrix.length;
        int numCols = Matrix[0].length;
        double[][] temp_Matrix = new double[numRows][numCols];
        for (int i = 0; i < numRows; i++) {
            for (int k = 0; k < numCols; k++) {
                temp_Matrix[i][k] = Matrix[i][k];
            }
        }
        for (int i = 0; i < temp.length; i++) {
            temp_Matrix[row_index][i] = temp[i] * factor;
        }
        return temp_Matrix;
    }
    public static void Print_Vector(double[] Vector) {
        for (int i = 0; i < Vector.length; i++) {
            System.out.print(Vector[i] + " ");
        }
        System.out.println("");
    }
    public static void Print_Vector(int[] Vector) {
        for (int i = 0; i < Vector.length; i++) {
            System.out.print(Vector[i] + " ");
        }
        System.out.println("");
    }
    public static String getSolutionType(double[][] RREF_Matrix) {
    int numRows = RREF_Matrix.length;
    int numCols = RREF_Matrix[0].length;
    int numVariables = numCols - 1; // Assuming it's an augmented matrix

    // Check for no solution
    for (int i = 0; i < numRows; i++) {
        boolean allZerosExceptLast = true;
        for (int j = 0; j < numVariables; j++) {
            if (RREF_Matrix[i][j] != 0.0) {
                allZerosExceptLast = false;
                break;
            }
        }
        if (allZerosExceptLast && RREF_Matrix[i][numVariables] != 0.0) {
            return "No Solution";
        }
    }

    // Check for free variables
    int pivotCount = 0;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numVariables; j++) {
            if (RREF_Matrix[i][j] == 1.0) {
                // Check if it's a leading 1
                boolean isLeadingOne = true;
                for (int k = 0; k < j; k++) {
                    if (RREF_Matrix[i][k] != 0.0) {
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
    
    private static double rounding(double num, int decimal){
        return Math.round(Math.pow(10, decimal)*num)/Math.pow(10, decimal);
    }

    public static double[][] RREF(double[][] Matrix, boolean debug) {
    int numRows = Matrix.length;
    int numCols = Matrix[0].length;
    double[][] temp_Matrix = new double[numRows][numCols];

    for (int i = 0; i < numRows; i++) {
        for (int k = 0; k < numCols; k++) {
            temp_Matrix[i][k] = Matrix[i][k];
        }
    }


    int pivotRow = 0;

    for (int j = 0; j < numCols && pivotRow < numRows; j++) {
        int maxRow = pivotRow;
        for (int i = pivotRow + 1; i < numRows; i++) {
            if (Math.abs(temp_Matrix[i][j]) > Math.abs(temp_Matrix[maxRow][j])) {
                maxRow = i;
            }
        }

        if (Math.abs(temp_Matrix[maxRow][j]) > epsilon) { 
            temp_Matrix = switch_rows(temp_Matrix, pivotRow, maxRow);
            if(debug){
                therefore();
                Print_Matrix(temp_Matrix);
            }

            double pivotValue = temp_Matrix[pivotRow][j];
            temp_Matrix = Constant_factor(temp_Matrix, pivotRow, 1/pivotValue);
            temp_Matrix = roundMatrix(temp_Matrix, epsilon);
            if(debug){
                therefore();
                Print_Matrix(temp_Matrix);
            }

            for(int i = pivotRow + 1; i < numRows; i++){
                double factor = temp_Matrix[i][j];
                temp_Matrix = replacement(temp_Matrix, i, pivotRow, -factor);
                temp_Matrix = roundMatrix(temp_Matrix, epsilon);
                if(debug){
                    therefore();
                    Print_Matrix(temp_Matrix);
                }
            }
            pivotRow++;
        }
    }
    
    for(int i = numRows - 1; i >=0; i--){
        int pivotCol = -1;
        for(int j = 0; j < numCols; j++){
            if(Math.abs(temp_Matrix[i][j] - 1.0) < epsilon){ 
                pivotCol = j;
                break;
            }
        }
        if(pivotCol != -1){
            for(int k = 0; k < i; k++){
                double factor = temp_Matrix[k][pivotCol];
                temp_Matrix = replacement(temp_Matrix, k, i, -factor);
                temp_Matrix = roundMatrix(temp_Matrix, epsilon);
                if(debug){
                    therefore();
                    Print_Matrix(temp_Matrix);
                }
            }
        }
    }
    
    temp_Matrix = roundMatrix(temp_Matrix, epsilon);
    if(debug){
        therefore();
    }
    return temp_Matrix;
}
    public static double[][] RREF(double[][] Matrix) {
    int numRows = Matrix.length;
    int numCols = Matrix[0].length;
    double[][] temp_Matrix = new double[numRows][numCols];

    for (int i = 0; i < numRows; i++) {
        for (int k = 0; k < numCols; k++) {
            temp_Matrix[i][k] = Matrix[i][k];
        }
    }

    
    int pivotRow = 0;

    for (int j = 0; j < numCols && pivotRow < numRows; j++) {
        int maxRow = pivotRow;
        for (int i = pivotRow + 1; i < numRows; i++) {
            if (Math.abs(temp_Matrix[i][j]) > Math.abs(temp_Matrix[maxRow][j])) {
                maxRow = i;
            }
        }

        if (Math.abs(temp_Matrix[maxRow][j]) > epsilon) { 
            temp_Matrix = switch_rows(temp_Matrix, pivotRow, maxRow);

            double pivotValue = temp_Matrix[pivotRow][j];
            temp_Matrix = Constant_factor(temp_Matrix, pivotRow, 1/pivotValue);
            temp_Matrix = roundMatrix(temp_Matrix, epsilon);

            for(int i = pivotRow + 1; i < numRows; i++){
                double factor = temp_Matrix[i][j];
                temp_Matrix = replacement(temp_Matrix, i, pivotRow, -factor);
                temp_Matrix = roundMatrix(temp_Matrix, epsilon);
            }
            pivotRow++;
        }
    }
    
    for(int i = numRows - 1; i >=0; i--){
        int pivotCol = -1;
        for(int j = 0; j < numCols; j++){
            if(Math.abs(temp_Matrix[i][j] - 1.0) < epsilon){ 
                pivotCol = j;
                break;
            }
        }
        if(pivotCol != -1){
            for(int k = 0; k < i; k++){
                double factor = temp_Matrix[k][pivotCol];
                temp_Matrix = replacement(temp_Matrix, k, i, -factor);
                temp_Matrix = roundMatrix(temp_Matrix, epsilon);
            }
        }
    }
    
    temp_Matrix = roundMatrix(temp_Matrix, epsilon);
    return temp_Matrix;
}
    private static double[][] roundMatrix(double[][] matrix, double epsilon) {
    int numRows = matrix.length;
    int numCols = matrix[0].length;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            if (Math.abs(matrix[i][j]) < epsilon) {
                matrix[i][j] = 0.0;
            }
        }
    }
    return matrix;
}
    public static void therefore() {
        System.out.println("->");
    }

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

    public static double[] matrixVectorMult(double[][] matrix, double[] vector){
        if(matrix[0].length != vector.length){
            throw new IllegalArgumentException(ColorText.errorFormat("Incompatible matrix and vector dimensions"));
        }

        double[] result = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            result[i] = dot_product(matrix[i], vector);
        }
        return result;
    }
    
    public static double dot_product(double[] vector1,double[] vector2, boolean debug){
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
    public static double nanoToSeconds(long nanoseconds){
        return nanoseconds/1e9;
    }
    public static double dot_product(double[] vector1,double[] vector2){
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
    public static double[][] generate_Matrix(int Rows, int Columns){
        double[][] Matrix = new double[Rows][Columns];
        Random random = new Random();
        for(int i = 0; i < Rows; i++){
            for(int k = 0; k < Columns; k++){
                Matrix[i][k] = random.nextInt(10);
            }
        }
        return Matrix;

    }
    public static double[][] Matrix_Mult(double[][] Matrix_A, double[][] Matrix_B, boolean debug){
        if(Matrix_A[0].length == Matrix_B.length){
            // proceed as usual type shit
            double[][] Result_Matrix = new double[Matrix_A.length][Matrix_B[0].length];
            for(int i = 0; i < Matrix_A.length; i++){
                for(int j = 0; j < Matrix_B[0].length; j++){
                    double[] Row = access_row(Matrix_A, i);
                    double[] Column = access_column(Matrix_B, j);
                    Result_Matrix[i][j] = dot_product(Row, Column, debug);
                }
            }
            return Result_Matrix;
        } else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input Matrices are invalid"));

        }
    }
    public static double[][] Matrix_Mult(double[][] Matrix_A, double[][] Matrix_B){
        if(Matrix_A[0].length == Matrix_B.length){
            // proceed as usual type shit
            double[][] Result_Matrix = new double[Matrix_A.length][Matrix_B[0].length];
            for(int i = 0; i < Matrix_A.length; i++){
                for(int j = 0; j < Matrix_B[0].length; j++){
                    double[] Row = access_row(Matrix_A, i);
                    double[] Column = access_column(Matrix_B, j);
                    Result_Matrix[i][j] = dot_product(Row, Column);
                }
            }
            return Result_Matrix;
        } else{
            throw new IllegalArgumentException(ColorText.errorFormat("Input Matrices are invalid"));

        }
    }
    public static double[] cross_product(double[] vector1, double[] vector2){
        if(vector1.length != 3 || vector2.length != 3)
            throw new IllegalArgumentException("");

        double[] result_vector = {
            (vector1[2]*vector2[3]-vector1[3]*vector2[2]),
            (vector1[3]*vector2[1]-vector1[1]*vector2[3]),
            (vector1[1]*vector2[2]-vector1[2]*vector2[1])
        };
        return result_vector;
    }
    public static double[][] matrixTransposition(double[][] matrix){
        double[][] tempMatrix = new double[matrix[0].length][matrix.length];
        double[] tempRow = new double[tempMatrix[0].length];
        for(int i = 0; i<matrix.length;i++){
            tempRow = matrix[i].clone();
            for(int k = 0; k < tempMatrix.length; k++){
                tempMatrix[k][i] = tempRow[k];
            }
        }
        return tempMatrix;
    }
    public static double[][] subMatrix(double[][] matrix, int column) {
        double[][] tempMatrix = new double[matrix.length - 1][matrix[0].length - 1];
        for (int i = 1; i < matrix.length; i++) {
            int subCol = 0;
            for (int j = 0; j < matrix[0].length; j++) {
                if (j == column) continue;
                tempMatrix[i - 1][subCol] = matrix[i][j];
                subCol++;
            }
        }
        return tempMatrix;
    }


    public static double determinant(double[][] matrix){
		double sum = 0;
		if(matrix.length != matrix[0].length){
			throw new IllegalArgumentException(ColorText.errorFormat("determinent input is not equal"));
		}
		int n = matrix.length;
		
		if(n==1){
			return matrix[0][0];
		}
		if(n == 2){
			return (matrix[0][0] * matrix[1][1]-
			matrix[0][1] * matrix[1][0]);
		}
        double sign;
        for(int i = 0; i < n; i++){
            if(i%2==0){
                sign = 1.0;
            }else{
                sign = -1.0;
            }
            sum += sign* matrix[0][i] * determinant(subMatrix(matrix, i));
        }
        return sum;
	}

    public static double[][] inverseMatrix(double[][] Matrix){
        if(Matrix.length != Matrix[0].length){
            throw new IllegalArgumentException(ColorText.errorFormat("Not a square matrix"));
        }
        // big booty bitches ;)
        if(determinant(Matrix) == 0){
            throw new IllegalArgumentException(ColorText.errorFormat("Determinent is zero"));
        }
        int n = Matrix.length;

        if(n == 2){
            double det = determinant(Matrix);
            double[][] inverse = {
                {Matrix[1][1]/det, -Matrix[0][1]/det},
                {-Matrix[1][0]/det, Matrix[0][0]/det}
            };
            return inverse;
        }
        int numRows = Matrix.length;
        int numCols = Matrix[0].length;
        double[][] temp_Matrix = new double[numRows][numCols];
        double[][] inverse_Matrix = new double[numRows][numCols];
        for(int i = 0; i < numRows; i++) {
            for (int k = 0; k < numCols; k++) {
                temp_Matrix[i][k] = Matrix[i][k];
            }
        }

        for(int i = 0; i < numRows;i++){
            inverse_Matrix[i][i] = 1;
        }


        
        int pivotRow = 0;

        for (int j = 0; j < numCols && pivotRow < numRows; j++) {
            int maxRow = pivotRow;
            for (int i = pivotRow + 1; i < numRows; i++) {
                if (Math.abs(temp_Matrix[i][j]) > Math.abs(temp_Matrix[maxRow][j])) {
                    maxRow = i;
                }
            }

            if (Math.abs(temp_Matrix[maxRow][j]) > epsilon) { 
                temp_Matrix = switch_rows(temp_Matrix, pivotRow, maxRow);
                inverse_Matrix = switch_rows(inverse_Matrix, pivotRow, maxRow);

                double pivotValue = temp_Matrix[pivotRow][j];
                temp_Matrix = Constant_factor(temp_Matrix, pivotRow, 1/pivotValue);
                inverse_Matrix = Constant_factor(inverse_Matrix, pivotRow, 1/pivotValue);
                temp_Matrix = roundMatrix(temp_Matrix, epsilon);
                inverse_Matrix = roundMatrix(inverse_Matrix, epsilon);

                for(int i = pivotRow + 1; i < numRows; i++){
                    double factor = temp_Matrix[i][j];
                    temp_Matrix = replacement(temp_Matrix, i, pivotRow, -factor);
                    temp_Matrix = roundMatrix(temp_Matrix, epsilon);
                    inverse_Matrix = replacement(inverse_Matrix, i, pivotRow, -factor);
                    inverse_Matrix = roundMatrix(inverse_Matrix, epsilon);
                }
                pivotRow++;
            }
        }
        
        for(int i = numRows - 1; i >=0; i--){
            int pivotCol = -1;
            for(int j = 0; j < numCols; j++){
                if(Math.abs(temp_Matrix[i][j] - 1.0) < epsilon){ 
                    pivotCol = j;
                    break;
                }
            }
            if(pivotCol != -1){
                for(int k = 0; k < i; k++){
                    double factor = temp_Matrix[k][pivotCol];
                    temp_Matrix = replacement(temp_Matrix, k, i, -factor);
                    temp_Matrix = roundMatrix(temp_Matrix, epsilon);
                    inverse_Matrix = replacement(inverse_Matrix, k, i, -factor);
                    inverse_Matrix = roundMatrix(inverse_Matrix, epsilon);
                }
            }
        }
        
        inverse_Matrix = roundMatrix(inverse_Matrix, epsilon);
        return inverse_Matrix;
    }



    public static double[][] inverseMatrix(double[][] Matrix,boolean stopwatch){
        long startTime = System.nanoTime();
        
        if(Matrix.length != Matrix[0].length){
            throw new IllegalArgumentException(ColorText.errorFormat("Not a square matrix"));
        }
        // big booty bitches ;)
        if(determinant(Matrix) == 0){
            throw new IllegalArgumentException(ColorText.errorFormat("Determinent is zero"));
        }
        int n = Matrix.length;

        if(n == 2){
            
        }
        int numRows = Matrix.length;
        int numCols = Matrix[0].length;
        double[][] temp_Matrix = new double[numRows][numCols];
        double[][] inverse_Matrix = new double[numRows][numCols];
        for(int i = 0; i < numRows; i++) {
            for (int k = 0; k < numCols; k++) {
                temp_Matrix[i][k] = Matrix[i][k];
            }
        }

        for(int i = 0; i < numRows;i++){
            inverse_Matrix[i][i] = 1;
        }


        int pivotRow = 0;

        for (int j = 0; j < numCols && pivotRow < numRows; j++) {
            int maxRow = pivotRow;
            for (int i = pivotRow + 1; i < numRows; i++) {
                if (Math.abs(temp_Matrix[i][j]) > Math.abs(temp_Matrix[maxRow][j])) {
                    maxRow = i;
                }
            }

            if (Math.abs(temp_Matrix[maxRow][j]) > epsilon) { 
                temp_Matrix = switch_rows(temp_Matrix, pivotRow, maxRow);
                inverse_Matrix = switch_rows(inverse_Matrix, pivotRow, maxRow);

                double pivotValue = temp_Matrix[pivotRow][j];
                temp_Matrix = Constant_factor(temp_Matrix, pivotRow, 1/pivotValue);
                inverse_Matrix = Constant_factor(inverse_Matrix, pivotRow, 1/pivotValue);
                temp_Matrix = roundMatrix(temp_Matrix, epsilon);
                inverse_Matrix = roundMatrix(inverse_Matrix, epsilon);

                for(int i = pivotRow + 1; i < numRows; i++){
                    double factor = temp_Matrix[i][j];
                    temp_Matrix = replacement(temp_Matrix, i, pivotRow, -factor);
                    temp_Matrix = roundMatrix(temp_Matrix, epsilon);
                    inverse_Matrix = replacement(inverse_Matrix, i, pivotRow, -factor);
                    inverse_Matrix = roundMatrix(inverse_Matrix, epsilon);
                }
                pivotRow++;
            }
        }
        
        for(int i = numRows - 1; i >=0; i--){
            int pivotCol = -1;
            for(int j = 0; j < numCols; j++){
                if(Math.abs(temp_Matrix[i][j] - 1.0) < epsilon){ 
                    pivotCol = j;
                    break;
                }
            }
            if(pivotCol != -1){
                for(int k = 0; k < i; k++){
                    double factor = temp_Matrix[k][pivotCol];
                    temp_Matrix = replacement(temp_Matrix, k, i, -factor);
                    temp_Matrix = roundMatrix(temp_Matrix, epsilon);
                    inverse_Matrix = replacement(inverse_Matrix, k, i, -factor);
                    inverse_Matrix = roundMatrix(inverse_Matrix, epsilon);
                }
            }
        }
        
        inverse_Matrix = roundMatrix(inverse_Matrix, epsilon);
        if(stopwatch){
            long estimatedTime = System.nanoTime() - startTime;
            System.out.println("elapsed time in seconds: "+nanoToSeconds(estimatedTime));
        }
        return inverse_Matrix;
    }

    }   