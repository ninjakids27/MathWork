public class Tensor {
    
    private int[] shape;
    private int dimensions;
    private double[] data;


    public Tensor(){
        this.shape = new int[]{3,3,3};
        this.dimensions = 3;
        this.data = new double[27];
    }
    public Tensor(int dimensions){
        this.shape = new int[dimensions];
        for(int i = 0; i < dimensions;i++){
            this.shape[i]=3;
        }
        this.dimensions = dimensions;
        this.data = new double[dimensions*3];
    }

}
