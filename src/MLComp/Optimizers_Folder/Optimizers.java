public class Optimizers {
    
    public static void gradientDescent(Neuron neuron,double[] gradient, double learningRate){
        // implement algorithm here
        double[] weights = neuron.getWeights();
        for(int i = 0; i < weights.length; i++){
            // Update rule: weight = weight - (learningRate * gradient)
            weights[i] -= learningRate * gradient[i];
        }
        neuron.setWeights(weights);
    }

    public static void sgd(Neuron neuron, double[] gradient, double learningRate) {
        double[] weights = neuron.getWeights();
        for (int i = 0; i < weights.length; i++) {
            // Update rule: weight = weight - (learningRate * gradient)
            weights[i] -= learningRate * gradient[i];
        }
        
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
    public static void sgdMomentum(Neuron neuron, double[] gradients, double learningRate) {
        double momentum = 0.9;
        double[] weights = neuron.getWeights();
        // Initialize the velocity array if it's the first pass
        if (velocity == null || velocity.length != weights.length) {
            velocity = new double[weights.length];
        }

        for (int i = 0; i < weights.length; i++) {
            // Update rule: v = (momentum * v) + (lr * gradient)
            velocity[i] = (momentum * velocity[i]) + (learningRate * gradients[i]);
            
            // Apply update: weight = weight - v
            weights[i] -= velocity[i];
        }
    }
    
    public static void adaGrad(Neuron neuron,double[] gradident,  double learningRate){
        
    }
    
    public static void adam(Neuron neuron,double[] gradident){
        
    }
}
