public class Optimizers {
    
    public static void gradientDescent(Neuron neuron,double[] gradident){
        neuron.setWeights(MatrixOps.vectorSum(MatrixOps.vectorMult(-1, gradident), neuron.getWeights()));
        neuron.setBias(neuron.getBias() + gradident[gradident.length - 1]);
    }

    public static void sgd(Neuron neuron, double[] gradient) {
    // 1. Get current weights and bias
    double[] weights = neuron.getWeights();
    double bias = neuron.getBias();
    
    // 2. Update the weights
    // We loop up to weights.length.
    // This assumes the gradient array contains [w1, w2, ... wn, bias]
    for (int i = 0; i < weights.length; i++) {
        weights[i] -= gradient[i];
    }
    
    // 3. Update the bias using the very last element of the gradient array
    // gradient.length should be weights.length + 1
    bias -= gradient[gradient.length - 1];
    
    // 4. Save the updated values back to the neuron
    neuron.setWeights(weights);
    neuron.setBias(bias);
}
    
    public static void sgdMomentum(Neuron neuron,double[] gradident){
        
    }
    
    public static void adaGrad(Neuron neuron,double[] gradident,  double learningRate){
        
    }
    
    public static void adam(Neuron neuron,double[] gradident){
        
    }
}
