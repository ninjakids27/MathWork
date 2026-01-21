package MathComp;


public class Tensor<T> {
    private final T[] data;
    private final int[] dimensions;
    private final int rank;
    // ??? 
    // this is NEVER going to get done - Programmer 2025
    public Tensor(T[] data, int[] dimensions){
        this.dimensions = dimensions;
        this.data = data;
        this.rank = dimensions.length;
    }

    @SuppressWarnings("unchecked")
    public Tensor(int[] dimensions){
        this.dimensions = dimensions;
        int mult =1;
        for(int dimension : dimensions){
            mult*=dimension;
        }
        this.data = (T[]) new Object[mult];
        this.rank = dimensions.length;
    }


   public T get(int[] coordinates){
        int index = 0;
        int base = 1;

        for (int i = coordinates.length - 1; i >= 0; i--) {
            index += coordinates[i] * base;
            
            base *= dimensions[i]; 
        }
        return data[index];
   }
   public int[] getDimensions(){
    return this.dimensions;
   }
   

}
