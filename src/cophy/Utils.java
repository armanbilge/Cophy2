/**
 * Utils.java
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

package cophy;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.MachineAccuracy;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public final class Utils {

    private Utils() {}

    public static final List<Node> getLineagesInHeightRange(Tree tree,
            double lower, double upper) {
        
        if (lower > upper)
            throw new IllegalArgumentException("Must be lower <= upper");
        
        List<Node> lineages = new ArrayList<Node>(tree.getLeafNodeCount());
        getLineagesInHeightRange(tree.getRoot(), lower, upper, lineages);

        return lineages;
        
    }
    
    private static final void getLineagesInHeightRange(Node node,
            double lower, double upper, List<Node> lineages) {
        
        if (node.isRoot() || node.getParent().getHeight() >= lower) {
            if (node.getHeight() < upper) lineages.add(node);
        } else {
            return;
        }
        
        for (Node child : node.getChildren())
            getLineagesInHeightRange(child, lower, upper, lineages);
        
    }
    
    public static final boolean lineageExistedAtHeight(Node node,
            double height) {
        
        return node.getHeight() <= height &&
                (node.isRoot() || node.getParent().getHeight() > height);
        
    }
    
    public static final List<Node> getLineagesAtHeight(Tree tree,
            double height) {
        return getLineagesAtHeight(tree, height, false);
    }
    
    public static final List<Node> getLineagesAtHeight(Tree tree, double height,
            boolean approachFromPresent) {
        
        List<Node> lineages = new ArrayList<Node>();
        if (approachFromPresent)
            getLineagesAtHeightApproachingFromPresent(tree.getRoot(),
                    height, lineages);
        else
            getLineagesAtHeight(tree.getRoot(), height, lineages);
        return lineages;
        
    }
    
    private static final void getLineagesAtHeight(Node node, double height,
            List<Node> lineages) {
                
        if (lineageExistedAtHeight(node, height)) {
            lineages.add(node);
        } else {
            getLineagesAtHeight(node.getLeft(), height, lineages);
            getLineagesAtHeight(node.getRight(), height, lineages);
        }
        
    }

    private static final void getLineagesAtHeightApproachingFromPresent(
            Node node, double height, List<Node> lineages) {
                
        if (node.getHeight() < height) {
            lineages.add(node);
        } else if (!node.isLeaf()) {
            getLineagesAtHeightApproachingFromPresent(node.getLeft(), height,
                    lineages);
            getLineagesAtHeightApproachingFromPresent(node.getRight(), height,
                    lineages);
        }
        
    }

    public static final int getLineageCountAtHeight(Tree tree,
            double height, boolean approachFromPresent) {
        
        if (approachFromPresent)
            return getLineageCountAtHeightApproachingFromPresent(tree.getRoot(),
                    height);
        else
            return getLineageCountAtHeight(tree.getRoot(), height);
        
    }
    
    private static final int getLineageCountAtHeight(Node node, double height) {
                
        if (lineageExistedAtHeight(node, height)) {
            return 1;
        } else {
            return getLineageCountAtHeight(node.getLeft(), height)
                + getLineageCountAtHeight(node.getRight(), height);
        }
        
    }

    private static final int getLineageCountAtHeightApproachingFromPresent(
            Node node, double height) {
                
        if (node.getHeight() < height &&
                (node.isRoot() || node.getParent().getHeight() > height)) {
            return 1;
        } else {
            return getLineageCountAtHeight(node.getLeft(), height)
                + getLineageCountAtHeight(node.getRight(), height);
        }
        
    }
    
    public static enum Relationship {
        SELF,
        DESCENDANT,
        ANCESTOR,
        SISTER,
        COUSIN;
    }
    
    public static final Relationship determineRelationship(Node self,
            Node relative) {
        
        if (self.equals(relative)) return Relationship.SELF;
        
        if (!self.isRoot() && !relative.isRoot()
                && self.getParent().equals(relative.getParent()))
            return Relationship.SISTER;
        
        Node intermediate = relative;
        while (!intermediate.isRoot()) {
            intermediate = intermediate.getParent();
            if (intermediate.equals(self)) return Relationship.DESCENDANT;
        }
        
        intermediate = self;
        while (!intermediate.isRoot()) {
            intermediate = intermediate.getParent();
            if (intermediate.equals(relative)) return Relationship.ANCESTOR;
        }
        
        return Relationship.COUSIN;
    }
    
    public static final boolean isCompleteCherry(Node parent, Node left,
            Node right) {
        
        return determineRelationship(left, right) == Relationship.SISTER
                && left.getParent().equals(parent);
        
    }
    
    public static final int sum(int[] values) {
        
        int sum = 0;
        for (int v : values) sum += v;
        return sum;
        
    }
    
    /**
     * Uses MachineAccuracy.same() to check equivalence. Note: this comparator
     * imposes orderings that are inconsistent with equals.
     */
    public static final class RelaxedComparator implements Comparator<Double> {
        @Override
        public int compare(Double d1, Double d2) {
            if (MachineAccuracy.same(d1, d2))
                return 0;
            else if (d1 > d2)
                return 1;
            else
                return -1;
        }
    }

    
}
