public class Optimizers {
    
    public static void gradientDescent(Neuron neuron,double[] gradident){
        neuron.setWeights(MatrixOps.vectorSum(MatrixOps.vectorMult(-1, gradident), neuron.getWeights()));
        neuron.setBias(neuron.getBias() - gradident[gradident.length - 1]);
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
    


    // You need a way to store the previous velocities between calls.
    // In a multi-layer network, you'll eventually want a Map or a way 
    // to track this per layer/weight-set.
    private static double[] velocity;

    /**
     * SGD with Momentum
     * @param weights The current weights to be updated
     * @param gradients The gradients calculated during backprop
     * @param lr The learning rate (e.g., 0.01)
     * @param momentum The momentum coefficient (usually 0.9)
     */
    public static void sgdMomentum(Neuron neuron, double[] gradients) {
        double lr = 0.01;
        double momentum = 0.9;
        double[] weights = neuron.getWeights();
        // Initialize the velocity array if it's the first pass
        if (velocity == null || velocity.length != weights.length) {
            velocity = new double[weights.length];
        }

        for (int i = 0; i < weights.length; i++) {
            // Update rule: v = (momentum * v) + (lr * gradient)
            velocity[i] = (momentum * velocity[i]) + (lr * gradients[i]);
            
            // Apply update: weight = weight - v
            weights[i] -= velocity[i];
        }
    }
    
    public static void adaGrad(Neuron neuron,double[] gradident,  double learningRate){
        
    }
    
    public static void adam(Neuron neuron,double[] gradident){
        
    }
}
