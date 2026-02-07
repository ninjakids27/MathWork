package MathComp;


public class Tensor<T>{
    protected final T[] data;
    private final int[] dimensions;
    private final int rank;
    // ??? 
    // this is NEVER going to get done - Programmer 2025
    public Tensor(T[] data, int[] dimensions){
        this.dimensions = dimensions;
        this.data = data;
        this.rank = dimensions.length;
    }
    // Copy constructor
    public Tensor(Tensor<T> other){
        this.dimensions = other.getDimensions();
        this.data = other.data;
        this.rank = other.rank;
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


   public void set(int[] coordinates, T value){
    int index = 0;
        int base = 1;

        for (int i = coordinates.length - 1; i >= 0; i--) {
            index += coordinates[i] * base;
            
            base *= dimensions[i]; 
        }
    data[index] = value;
   }
   public int[] getDimensions(){
    return this.dimensions;
   }
   

}
