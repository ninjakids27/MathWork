package MathComp;

public class Matrix<T> extends Tensor<T>{
    

    public Matrix(T[] data,int rows, int columns){
        super(data,new int[] {rows,columns});
    }

    @Override
    public String toString(){
        int num_rows = super.getDimensions()[0];
        int num_columns = super.getDimensions()[1];
        String a = "";
        for (int i = 0; i < num_rows; i++) {
            for (int k = 0; k < num_columns; k++) {
                a+=super.get(new int[] {i,k})+" ";
            }
            a+="\n";
        }
        return a;
    }
}
