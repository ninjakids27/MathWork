package MLComp.Optimizers_Folder;

import MLComp.Neuron;
import MathComp.MatrixOps;
import Runner.ColorText;

public class Optimizers {

    public static void sgd(Neuron neuron, double[] gradient, double learningRate) {
        if(neuron.getWeightsLength() + 1 != gradient.length){
            MatrixOps.printVector(gradient);
            double[] weights = neuron.getWeights();
            MatrixOps.printVector(weights);
            System.out.println(neuron.getBias());
            throw new IllegalArgumentException(ColorText.errorFormat("Gradient length does not match weights + bias length."));
        }
        double[] weights = neuron.getWeights();
        for (int i = 0; i < weights.length; i++) {
            // Update rule: weight = weight - (learningRate * gradient)
            weights[i] -= learningRate * gradient[i];
        }
        neuron.setWeights(weights);

        // Update bias (gradient for bias is at the end of the gradient array)
        double bias = neuron.getBias();
        bias -= learningRate * gradient[gradient.length - 1];
        neuron.setBias(bias);
    }

    // You need a way to store the previous velocities between calls.
    // In a multi-layer network, you'll eventually want a Map or a way
    // to track this per layer/weight-set.
    private static double[] velocity;

    /**
     * SGD with Momentum
     * 
     * @param weights   The current weights to be updated
     * @param gradients The gradients calculated during backprop
     * @param lr        The learning rate (e.g., 0.01)
     * @param momentum  The momentum coefficient (usually 0.9)
     */
    public static void sgdMomentum(Neuron neuron, double[] gradients, double learningRate) {
        double momentum = 0.9;
        double[] weights = neuron.getWeights();
        double bias = neuron.getBias();

        // Initialize the velocity array if it's the first pass
        if (velocity == null || velocity.length != weights.length + 1) {
            velocity = new double[weights.length + 1];
        }

        for (int i = 0; i < weights.length; i++) {
            // Update rule: v = (momentum * v) + (lr * gradient)
            velocity[i] = (momentum * velocity[i]) + (learningRate * gradients[i]);

            // Apply update: weight = weight - v
            weights[i] -= velocity[i];
        }
        neuron.setWeights(weights);

        // Update bias
        int biasIndex = weights.length;
        velocity[biasIndex] = (momentum * velocity[biasIndex]) + (learningRate * gradients[biasIndex]);
        bias -= velocity[biasIndex];
        neuron.setBias(bias);
    }

    public static void adaGrad(Neuron neuron, double[] gradient, double learningRate) {
        double epsilon = 1e-8;
        double[] weights = neuron.getWeights();
        double bias = neuron.getBias();

        // Initialize sum of squared gradients (V) if first pass or size mismatch
        if (neuron.getV() == null || neuron.getV().length != weights.length) {
            neuron.setV(new double[weights.length]);
            neuron.setVBias(0);
        }

        double[] sumSquaredGradients = neuron.getV();

        // Update weights
        for (int i = 0; i < weights.length; i++) {
            // Accumulate squared gradients
            sumSquaredGradients[i] += gradient[i] * gradient[i];

            // Update weight: w = w - (lr / sqrt(G + eps)) * g
            weights[i] -= (learningRate / (Math.sqrt(sumSquaredGradients[i]) + epsilon)) * gradient[i];
        }
        neuron.setWeights(weights);

        // Update bias
        // The bias gradient is stored at the end of the gradient array
        double biasGradient = gradient[gradient.length - 1];

        double vBias = neuron.getVBias();
        vBias += biasGradient * biasGradient;
        neuron.setVBias(vBias);

        bias -= (learningRate / (Math.sqrt(vBias) + epsilon)) * biasGradient;
        neuron.setBias(bias);
    }

    /**
     * Adam Optimizer
     * 
     * @param neuron       The neuron whose weights will be updated
     * @param gradient     The gradients calculated during backprop
     * @param learningRate The learning rate (e.g., 0.001)
     */
    public static void adam(Neuron neuron, double[] gradient, double learningRate) {
        double beta1 = 0.9;
        double beta2 = 0.999;
        double epsilon = 1e-8;

        double[] weights = neuron.getWeights();
        double bias = neuron.getBias();

        // Initialize moment arrays if first pass or size mismatch
        if (neuron.getM() == null || neuron.getM().length != weights.length) {
            neuron.setM(new double[weights.length]);
            neuron.setV(new double[weights.length]);
            neuron.setMBias(0);
            neuron.setVBias(0);
            neuron.setT(0);
        }

        neuron.setT(neuron.getT() + 1); // Increment timestep
        int t = neuron.getT();

        // --- Update Weights ---
        for (int i = 0; i < weights.length; i++) {
            // Update biased first moment estimate
            neuron.getM()[i] = beta1 * neuron.getM()[i] + (1 - beta1) * gradient[i];

            // Update biased second raw moment estimate
            neuron.getV()[i] = beta2 * neuron.getV()[i] + (1 - beta2) * gradient[i] * gradient[i];

            // Compute bias-corrected first moment estimate
            double mHat = neuron.getM()[i] / (1 - Math.pow(beta1, t));

            // Compute bias-corrected second raw moment estimate
            double vHat = neuron.getV()[i] / (1 - Math.pow(beta2, t));

            // Update weights
            weights[i] -= learningRate * mHat / (Math.sqrt(vHat) + epsilon);
        }
        neuron.setWeights(weights);

        // --- Update Bias ---
        // The bias gradient is stored at the end of the gradient array
        double biasGradient = gradient[gradient.length - 1];

        // Update biased first moment estimate for bias
        neuron.setMBias(beta1 * neuron.getMBias() + (1 - beta1) * biasGradient);

        // Update biased second raw moment estimate for bias
        neuron.setVBias(beta2 * neuron.getVBias() + (1 - beta2) * biasGradient * biasGradient);

        // Compute bias-corrected first moment estimate for bias
        double mHatBias = neuron.getMBias() / (1 - Math.pow(beta1, t));

        // Compute bias-corrected second raw moment estimate for bias
        double vHatBias = neuron.getVBias() / (1 - Math.pow(beta2, t));

        // Update bias
        neuron.setBias(bias - learningRate * mHatBias / (Math.sqrt(vHatBias) + epsilon));
    }

}