/**
 * DHSLSim.java
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

package cophy.simulation;

import java.util.List;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import cophy.Util;

/**
 * 
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class DHSLSimulator {

    public static final String RECONCILIATION = "reconciliation";
    private static final String EXTINCT = "extinct";
    
    protected Tree host;
    protected double lambda;
    protected double tau;
    protected double mu;
    protected double origin;
    protected boolean complete;
    
    public DHSLSimulator(Tree host, double lambda, double tau, double mu,
            double origin, boolean complete) {
        this.host = host;
        this.lambda = lambda;
        this.tau = tau;
        this.mu = mu;
        this.origin = origin;
        this.complete = complete;
    }
    
    public Tree simulateTree() {
        return simulateTree(0.0);
    }
    
    public Tree simulateTree(double until) {
        return new Tree(simulateSubtree(host.getRoot(), origin, until));
    }
    
    public Tree resumeSimulation(Tree tree, double until) {
        
        for (final Node node : tree.getExternalNodes()) {
            
            if (!(Boolean) node.getMetaData(EXTINCT)) {
                
                final Node replacement = simulateSubtree(node,
                        (Node) node.getMetaData(RECONCILIATION),
                        node.getHeight(), until);
                
                if (!complete && replacement == null) {
                    if (node.isRoot()) {
                        return null;
                    } else {
                        final Node parent = node.getParent();
                        final Node sister = Util.isLeft(node) ?
                                parent.getRight() : parent.getLeft();
                        if (parent.isRoot()) {
                            tree.setRoot(sister);
                        } else {
                            final Node parentParent = parent.getParent();
                            final boolean parentIsLeft = Util.isLeft(parent);
                            if (parentIsLeft)
                                parentParent.setLeft(sister);
                            else
                                parentParent.setRight(sister);
                        }
                    }
                }
            }
        }
        
        return tree;
        
    }
    
    protected Node simulateSubtree(Node host, double height, double until) {
        return simulateSubtree(new Node(), host, height, until);
    }
    
    protected Node simulateSubtree(Node node, Node host, double height,
            double until) {
        
        node.setMetaData(RECONCILIATION, host);
        
        final Node left;
        final Node right;
        
        final int event;
        if (host.isRoot()) { // No host-switching at root
            event = Util.nextWeightedInteger(lambda, mu);
            height -= Util.nextPoissonTime(lambda, mu);
        } else {
            event = Util.nextWeightedInteger(lambda, mu, tau);
            height -= Util.nextPoissonTime(lambda, mu, tau);
        }
        
        final Node leftHost;
        final Node rightHost;
        final double hostHeight = host.getHeight();
        if (hostHeight > height) {
            height = hostHeight;
            
            if (host.isLeaf()) {
                node.setHeight(height);
                return node;
            }
            
            // Cospeciation event
            leftHost = host.getLeft();
            rightHost = host.getRight();
        } else if (until > height) {
            node.setHeight(until);
            return node;
        } else {
            switch(event) {
            case 0: // Duplication event
                leftHost = host;
                rightHost = host;
                break;
            case 1: // Loss event
                if (!complete) return null;
                node.setMetaData(EXTINCT, true);
                leftHost = null;
                rightHost = null;
                break;
            case 2: // Host-switch event
                final Node host1 = host;
                final List<Node> potentialHosts =
                        Util.getLineagesAtHeight(this.host, height);
                potentialHosts.remove(host);
                final Node host2 = Util.getRandomElement(potentialHosts);
                if (Randomizer.nextBoolean()) {
                    leftHost = host1;
                    rightHost = host2;
                } else {
                    leftHost = host1;
                    rightHost = host2;
                }
            default:
                throw new RuntimeException();
            }
        }
        
        if (leftHost == null) {
            left = null;
            right = null;
        } else {
            left = simulateSubtree(leftHost, height, until);
            right = simulateSubtree(rightHost, height, until);
        }
        
        node.setHeight(height);
        
        if (event != 1) {
            if (left != null && right != null) { // Both lineages survived
                node.setLeft(left);
                node.setRight(right);
            } else if (left == right) {
                // Both are null, hence this entire lineage was lost
                return null;
            } else { // Just one child lineage is null/lost
                return left != null ? left : right;
            }
        }
        
        return node;
    }
    
	public static Tree simulateTree(Tree host, double lambda, double tau,
	        double mu, double origin, boolean complete) {
		return new DHSLSimulator(host, lambda, tau, mu, origin, complete)
		    .simulateTree();
	}
		
}
