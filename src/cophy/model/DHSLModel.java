/**
 * DHSLModel.java
 * 
 * Cophy: Cophylogenetics for BEAST 2
 * 
 * Copyright (C) 2014 Arman D. Bilge <armanbilge@gmail.com>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package cophy.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.NavigableMap;
import java.util.Queue;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import cophy.Particle.TreeParticle;
import cophy.Reconciliation;
import cophy.Util;
import cophy.simulation.DHSLSimulator;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class DHSLModel extends EmbeddedTreeDistribution {

    public static final String RECONSTRUCTED = "reconstructed";
    
    public Input<Integer> particleCountInput =
            new Input<Integer>("particleCount",
            "The number of particles to simulate.",
            512);
    
    protected int particleCount;
    protected List<TreeParticle> particles;
    protected DHSLSimulator simulator; 
    
    public void initAndValidate() {
        super.initAndValidate();
        particleCount = particleCountInput.get();
        particles = new ArrayList<TreeParticle>(particleCount);
        IntStream.range(0, particleCount).forEach(i -> particles.add(null));
    }
    
    @Override
    protected double calculateLogDensity() {
        
        if (!isValid()) return Double.NEGATIVE_INFINITY;
        
        final Tree embeddedTree = embeddedTreeInput.get();
        final Reconciliation reconciliation = reconciliationInput.get();
        
        final double rate = rateParameterInput.get().getValue();
        simulator = new DHSLSimulator(hostTreeInput.get(),
                duplicationRateParameterInput.get().getValue() * rate,
                hostSwitchRateParameterInput.get().getValue() * rate,
                lossRateParameterInput.get().getValue() * rate,
                originHeightParameterInput.get().getValue(),
                true);
        
        generateParticles();
        
        particles.stream().forEach(particle -> {
            
            final Tree completeTree = particle.getValue();
            
            if (embeddedTree.getLeafNodeCount()
                    != Util.getLineageCountAtHeight(completeTree, 0.0, false)) {
                particle.muliWeight(0.0);
                return;
            }
            
            embeddedTree.getExternalNodes().stream().forEach(leaf -> {
                
                final Node host = reconciliation.getHost(leaf);
                
                final Node completeParent =
                        particle.getNodeMap().get(leaf.getParent());
                List<Node> potentialCompletes;
                
                final Node head;
                if (Util.isLeft(leaf))
                    head = completeParent.getLeft();
                else
                    head = completeParent.getRight();

                if (head.isLeaf())
                    potentialCompletes = Arrays.asList(head);
                else
                    potentialCompletes = head.getAllLeafNodes();
                
                potentialCompletes = potentialCompletes.stream()
                        .filter(n -> n.getMetaData(DHSLSimulator.RECONCILIATION)
                        .equals(host)).collect(Collectors.toList());
 
                final int size = potentialCompletes.size();
                if (size == 0) {
                    particle.muliWeight(0.0);
                    return;
                }
                
            });
            
        });
        
        final double p = particles.stream()
                .filter(particle -> particle.getWeight() > 0).count()
                / particleCount;
        
        return Math.log(p);
    }
    
    protected boolean isValid() {
        
        final Tree embedded = embeddedTreeInput.get();
        final Reconciliation reconciliation = reconciliationInput.get();
        
        for (int i = 0; i < embedded.getNodeCount(); ++i) {
            final Node node = embedded.getNode(i);
            final Node host = reconciliation.getHost(node);
            if (!Util.lineageExistedAtHeight(host, node.getHeight()))
                return false;
        }
        
        return true;
    }
        
    protected void generateParticles() {
        
        final Tree embedded = embeddedTreeInput.get();
        final Reconciliation reconciliation = reconciliationInput.get();
                
        final NavigableMap<Double,Node> heights2Nodes =
                new TreeMap<Double,Node>();
        for (Node node : embedded.getInternalNodes())
            heights2Nodes.put(node.getHeight(), node);
        final Queue<Double> speciationsQueue =
                new LinkedList<Double>(heights2Nodes.descendingKeySet());
                
        final double until = speciationsQueue.poll();
        IntStream.range(0, particleCount).forEach(i -> {
            final TreeParticle particle =
                    new TreeParticle(simulator.simulateTree(until));
            particles.set(i, particle);
            final Tree tree = particle.getValue();
            final Node reconstructed = heights2Nodes.get(until);
            assert(reconstructed.isRoot());
            final Node host = reconciliation.getHost(reconstructed);
            
            final List<Node> potentialCompletes =
                    tree.getExternalNodes().stream()
                    .filter(n -> n.getMetaData(DHSLSimulator.RECONCILIATION)
                            .equals(host)).collect(Collectors.toList());
            
            final int size = potentialCompletes.size();
            if (size == 0) {
                particle.muliWeight(0.0);
                return;
            }
            
            final Node complete = Util.getRandomElement(potentialCompletes);
            
            simulateSpeciation(particle, reconstructed, complete, host);
            particle.muliWeight(size);
        });

        Util.resample(particles);
        
        while (!speciationsQueue.isEmpty()) {
            
            final double until2 = speciationsQueue.poll();
            particles.stream().forEach(particle -> {
                Tree tree = particle.getValue();
                simulator.resumeSimulation(tree, until2);
                final Node reconstructed = heights2Nodes.get(until2);
                final Node host = reconciliation.getHost(reconstructed);
                
                final Node completeParent =
                        particle.getNodeMap().get(reconstructed.getParent());
                List<Node> potentialCompletes;
                
                final Node head;
                if (Util.isLeft(reconstructed))
                    head = completeParent.getLeft();
                else
                    head = completeParent.getRight();

                if (head.isLeaf())
                    potentialCompletes = Arrays.asList(head);
                else
                    potentialCompletes = head.getAllLeafNodes();
                
                potentialCompletes = potentialCompletes.stream()
                        .filter(n -> n.getMetaData(DHSLSimulator.RECONCILIATION)
                        .equals(host)).collect(Collectors.toList());
 
                final int size = potentialCompletes.size();
                if (size == 0) {
                    particle.muliWeight(0.0);
                    return;
                }
                
                final Node complete = Util.getRandomElement(potentialCompletes);
                
                simulateSpeciation(particle, reconstructed, complete, host);
                particle.muliWeight(size);

            });
            
            Util.resample(particles);
            
        }
        
        particles.stream().forEach(particle -> {
            final Tree tree = particle.getValue();
            simulator.resumeSimulation(tree, 0.0);
        });
        
    }
    
    protected void simulateSpeciation(final TreeParticle particle,
            final Node reconstructed, final Node complete, final Node host) {
        
        final double rate = rateParameterInput.get().getValue();
        final double lambda =
                duplicationRateParameterInput.get().getValue() * rate;
        final double tau =
                hostSwitchRateParameterInput.get().getValue() * rate;
        final double lambdaPtau = lambda + tau;
        
        final double height = reconstructed.getHeight();
        
        final Tree hostTree = hostTreeInput.get();
        
        particle.getNodeMap().put(reconstructed, complete);
        
        final Node child1 = new Node();
        final Node child2 = new Node();
        child1.setMetaData(DHSLSimulator.RECONCILIATION, host);
        child1.setHeight(height);
        child2.setHeight(height);
        
        if (Randomizer.nextDouble() < lambda / lambdaPtau) { // Duplication
            child2.setMetaData(DHSLSimulator.RECONCILIATION, host);
        } else { // Host-switch
            List<Node> potentialHosts = Util.getLineagesAtHeight(hostTree,
                    height);
            potentialHosts.remove(host);
            final int r = Randomizer.nextInt(potentialHosts.size());
            child2.setMetaData(DHSLSimulator.RECONCILIATION, potentialHosts.get(r));
        }
        
        if (Randomizer.nextBoolean()) {
            complete.setLeft(child1);
            complete.setRight(child2);
        } else {
            complete.setLeft(child2);
            complete.setRight(child1);
        }
        
        particle.muliWeight(lambdaPtau);
        
    }
    
}
