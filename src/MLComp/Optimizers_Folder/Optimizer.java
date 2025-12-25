// this is to make things look cleaner, I don't know if the optimizers can handle it
// since they can vary a lot so I don't know if this type of stuff is maintainable for
// optimizers, but I do know it's good for ActivaionFunctions.
@FunctionalInterface
public interface Optimizer {
    void update(Neuron Neuron,double[] gradident);
    
}