package cophy;

public class Particle<T> {

    T value;
    double weight;
    
    public Particle() {}
    
    public Particle(T value) {
        this(value, 1.0);
    }
    
    public Particle(T value, double weight) {
        this.value = value;
        this.weight = weight;
    }
    
    public T getValue() {
        return value;
    }
    
    public void setValue(T value) {
        this.value = value;
    }
    
    public double getWeight() {
        return weight;
    }
    
    public void setWeight(double weight) {
        this.weight = weight;
    }
    
}
