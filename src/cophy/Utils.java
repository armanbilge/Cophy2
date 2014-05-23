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
import java.util.List;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public final class Utils {

    public Utils() {}

    public static final List<Node> getLineagesInHeightRange(Tree tree,
            double lower, double upper) {
        
        if (lower > upper)
            throw new IllegalArgumentException("Must be lower <= upper");
        
        List<Node> lineages = new ArrayList<Node>(tree.getLeafNodeCount());
        getLineagesInHeightRange(tree.getRoot(), lower, upper,
                lineages);
        return lineages;
        
    }
    
    private static final void getLineagesInHeightRange(Node node,
            double lower, double upper, List<Node> lineages) {
        
        double nodeHeight = node.getHeight();
        if (nodeHeight >= lower) {
            if (nodeHeight < upper) lineages.add(node);
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
    
}
