package cophy.model;

import java.util.LinkedList;
import java.util.NavigableMap;
import java.util.Queue;
import java.util.TreeMap;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import cophy.Particle;
import cophy.Util;
import cophy.sim.DHSLSim;

public class DHSLModel extends EmbeddedTreeDistribution {

    public static final String RECONSTRUCTED = "reconstructed";
    
    public Input<Integer> particleCountInput =
            new Input<Integer>("particleCount",
            "The number of particles to simulate.",
            100);
    
    protected Particle<Tree>[] particles;
    protected DHSLSim simulator; 
    
    @SuppressWarnings("unchecked")
    public void initAndValidate() {
        super.initAndValidate();
        final int particleCount = particleCountInput.get();
        particles = new Particle[particleCount];
    }
    
    @Override
    protected double calculateLogDensity() {
        
        simulator = new DHSLSim(hostTreeInput.get(),
                duplicationRateParameterInput.get().getValue(),
                hostSwitchRateParameterInput.get().getValue(),
                lossRateParameterInput.get().getValue(),
                originHeightParameterInput.get().getValue(),
                true);
        
        Tree embedded = embeddedTreeInput.get();
        NavigableMap<Double,Node> speciations = new TreeMap<Double,Node>();
        for (Node node : embedded.getInternalNodes())
            speciations.put(node.getHeight(), node);
        Queue<Double> speciationsQueue =
                new LinkedList<Double>(speciations.descendingKeySet());
        
        double until = speciationsQueue.poll();
        for (int i = 0; i < particles.length; ++i) {
            particles[i] = new Particle<Tree>(simulator.simulateTree(until));
            
        }

        particles = Util.resample(particles);
        
        while (!speciationsQueue.isEmpty()) {
            until = speciationsQueue.poll();
            for (Particle<Tree> p : particles) {
                Tree tree = p.getValue();
                tree = simulator.resumeSimulation(tree, until);
                p.setValue(tree);
                if (tree == null) {
                    
                }
            }
        }
        
        return 0;
    }
    
    @Override
    protected double getOverallRate(Node embedded, Node host) {
        // TODO Auto-generated method stub
        return 0;
    }

}
