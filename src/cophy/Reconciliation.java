/**
 * Reconciliation.java
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

import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class Reconciliation extends IntegerParameter {

    protected Tree embeddedTree;
    protected Tree hostTree;
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        setLower(0);
        setUpper(hostTree.getNodeCount() - 1);
        setDimension(embeddedTree.getNodeCount());
    }
    
    public Node getHost(Node embeddedNode) {
        
        if (embeddedNode.getTree() != embeddedTree)
            throw new IllegalArgumentException("embeddedNode must belong to "
                    + "embeddedTree");

        return hostTree.getNode(getValue(embeddedNode.getNr()));
        
    }
    
    public void setHost(Node embeddedNode, Node hostNode) {
        
        if (embeddedNode.getTree() != embeddedTree)
            throw new IllegalArgumentException("embeddedNode must belong to "
                    + "embeddedTree");
        
        if (hostNode.getTree() != hostTree)
            throw new IllegalArgumentException("hostNode must belong to "
                    + "hostTree");

        setValue(embeddedNode.getNr(), hostNode.getNr());
        
    }
    
    public boolean notarizeTrees(Tree embeddedTree, Tree hostTree) {
        
        return embeddedTree == this.embeddedTree &&
                hostTree == this.hostTree;
    
    }
    
}
