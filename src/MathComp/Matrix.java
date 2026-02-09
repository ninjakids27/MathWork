package MathComp;

public class Matrix<T> extends Tensor<T>{
    private int[] coordinates = new int[] {0, 0};
    
    public T get(int row, int column) {
        coordinates[0] = row;
        coordinates[1] = column;
        return super.get(coordinates);
    }

    public void set(int row, int column, T value){
        coordinates[0] = row;
        coordinates[1] = column;
        super.set(coordinates, value);
    }

    public Matrix(T[] data,int rows, int columns){
        super(data,new int[] {rows,columns});
    }

    public Matrix(int rows, int columns,T[] data){
        super(data,new int[] {rows,columns});
    }

    public Matrix(Matrix<T> other){
        super(other.data,new int[] {other.getRowNum(),other.getColumnNum()});
    }


    public Matrix(int rows, int columns){
        super(new int[] {rows,columns});
    }

    public int getRowNum(){
        return super.getDimensions()[0];
    }
    
    
    public int getColumnNum(){
        return super.getDimensions()[1];
    }
    /**
     * Accesses a specific row of the given matrix.
     * @param Matrix the input matrix
     * @param row_index the index of the row to access (0-based)
     * @return the specified row as a 1D array
     */
    
    
    /**
     * Switches two rows in the given matrix and returns a new matrix with the rows swapped.
     * @param Matrix the input matrix
     * @param row1Index the index of the first row to swap (0-based)
     * @param row2Index the index of the second row to swap (0-based)
     * @return a new matrix with the specified rows swapped
     */
    public void switchRows(int row1Index, int row2Index) {
        // Clone the matrix
        T[] row1 = this.accessRow(row1Index);
        T[] row2 = this.accessRow(row2Index);
        for (int i = 0; i < row1.length; i++) {
            this.set(new int[] {row1Index,i}, row1[i]);
            this.set(new int[] {row2Index,i}, row2[i]);
        }
    }
    
    /**
     * Sets a specific row of the given matrix to the provided values.
     * @param rowNum the index of the row to set (0-based)
     * @param row the values to set in the row
     */
    public void setRow(int rowNum, T[] row){
        int[] coordinates = {rowNum,0};
        for(int i = 0; i < row.length; i++){
            coordinates[1] = i;
            this.set(coordinates, row[i]);
        }
    }
    /**
     * Accesses a specific column of the given matrix.
     * @param rowIndex the index of the column to access (0-based)
     * @return
     */
    @SuppressWarnings("unchecked")
    public T[] accessRow(int rowIndex) {
        int num_cols = this.getDimensions()[1];
        T[] temp = (T[]) new Object[num_cols];
        int[] coordinates =  {rowIndex,0};
        for (int i = 0; i < num_cols; i++) {
            temp[i] = this.get(coordinates);
            coordinates[1]++;
        }
        return temp;
    }


    @SuppressWarnings("unchecked")
    public T[] accessColumn(int columnIndex) {
        int num_rows = this.getDimensions()[0];
        T[] temp = (T[]) new Object[num_rows];
        int[] coordinates =  {0,columnIndex};
        for (int i = 0; i < num_rows; i++) {
            temp[i] = this.get(coordinates);
            coordinates[0]++;
        }
        return temp;
    }

     



    @Override
    public String toString(){
        int numRows = super.getDimensions()[0];
        int numColumns = super.getDimensions()[1];
        String a = "";
        for (int i = 0; i < numRows; i++) {
            for (int k = 0; k < numColumns; k++) {
                a+=super.get(new int[] {i,k})+" ";
            }
            a+="\n";
        }
        return a;
    }


}
