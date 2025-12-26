public class Optimizers {
    
    public static void gradientDescent(Neuron neuron,double[] gradient, double learningRate){
        double[] weights = neuron.getWeights();
        for(int i = 0; i < weights.length; i++){
            // Update rule: weight = weight - (learningRate * gradient)
            weights[i] -= learningRate * gradient[i];
        }
        // neuron.setWeights(weights);
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
    
    
    
    /**
     * Adam Optimizer
     * @param neuron The neuron whose weights will be updated
     * @param gradient The gradients calculated during backprop
     * @param learningRate The learning rate (e.g., 0.001)
     */
    public static void adam(Neuron neuron, double[] gradient, double learningRate){
        double beta1 = 0.9;    // Exponential decay rate for first moment
        double beta2 = 0.999;  // Exponential decay rate for second moment
        double epsilon = 1e-8; // Small constant to prevent division by zero
        
        double[] weights = neuron.getWeights();
        
        // Initialize moment arrays if first pass or size mismatch
        if (neuron.getM() == null || neuron.getM().length != weights.length) {
            neuron.setM(new double[weights.length]);
            neuron.setV(new double[weights.length]);
            neuron.setT(0);
        }
        
        neuron.setT(neuron.getT() + 1); // Increment timestep
        
        for (int i = 0; i < weights.length; i++) {
            // Update biased first moment estimate
            neuron.getM()[i] = beta1 * neuron.getM()[i] + (1 - beta1) * gradient[i];
            
            // Update biased second raw moment estimate
            neuron.getV()[i] = beta2 * neuron.getV()[i] + (1 - beta2) * gradient[i] * gradient[i];
            
            // Compute bias-corrected first moment estimate
            double mHat = neuron.getM()[i] / (1 - Math.pow(beta1, neuron.getT()));
            
            // Compute bias-corrected second raw moment estimate
            double vHat = neuron.getV()[i] / (1 - Math.pow(beta2, neuron.getT()));
            
            // Update weights
            weights[i] -= learningRate * mHat / (Math.sqrt(vHat) + epsilon);
        }
        neuron.setWeights(weights);

    }
}
