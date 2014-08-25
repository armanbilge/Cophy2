package cophy.sim;

import java.util.List;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import cophy.Util;

public class DHSLSim {

    private Tree host;
    private double lambda;
    private double tau;
    private double mu;
    private double origin;
    private boolean complete;
    
    private static final String RECONCILIATION = "reconciliation";
    
    public DHSLSim(Tree host, double lambda, double tau, double mu,
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
            
            final Node replacement = simulateSubtree(node,
                    (Node) node.getMetaData(RECONCILIATION), node.getHeight(),
                    until);
            
            if (replacement == null) {
                if (node.isRoot()) {
                    return null;
                } else {
                    final Node parent = node.getParent();
                    final Node sister = parent.getLeft() == node ?
                            parent.getRight() : parent.getLeft();
                    if (parent.isRoot()) {
                        tree.setRoot(sister);
                    } else {
                        final Node parentParent = parent.getParent();
                        final boolean parentIsLeft =
                                parentParent.getLeft() == parent;
                        if (parentIsLeft)
                            parentParent.setLeft(sister);
                        else
                            parentParent.setRight(sister);
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
                leftHost = null;
                rightHost = null;
                break;
            case 2: // Host-switch event
                leftHost = null;
                final List<Node> potentialHosts =
                        Util.getLineagesAtHeight(this.host, height);
                potentialHosts.remove(host);
                final int r = Randomizer.nextInt(potentialHosts.size());
                rightHost = potentialHosts.get(r);
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
		return new DHSLSim(host, lambda, tau, mu, origin, complete)
		    .simulateTree();
	}
		
}
