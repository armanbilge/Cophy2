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

import java.io.PrintStream;

import beast.core.StateNode;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class Reconciliation extends StateNode {

    protected Tree embeddedTree;
    protected Tree hostTree;
    protected int[] map;
    protected int[] storedMap;
    
    /**
     * 
     */
    public Reconciliation() {}

    /* (non-Javadoc)
     * @see beast.core.Loggable#init(java.io.PrintStream)
     */
    @Override
    public void init(PrintStream out) throws Exception {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see beast.core.Loggable#log(int, java.io.PrintStream)
     */
    @Override
    public void log(int nSample, PrintStream out) {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see beast.core.Loggable#close(java.io.PrintStream)
     */
    @Override
    public void close(PrintStream out) {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see beast.core.Function#getDimension()
     */
    @Override
    public int getDimension() {
        // TODO Auto-generated method stub
        return 0;
    }

    /* (non-Javadoc)
     * @see beast.core.Function#getArrayValue()
     */
    @Override
    public double getArrayValue() {
        // TODO Auto-generated method stub
        return 0;
    }

    /* (non-Javadoc)
     * @see beast.core.Function#getArrayValue(int)
     */
    @Override
    public double getArrayValue(int iDim) {
        // TODO Auto-generated method stub
        return 0;
    }

    /* (non-Javadoc)
     * @see beast.core.StateNode#setEverythingDirty(boolean)
     */
    @Override
    public void setEverythingDirty(boolean isDirty) {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see beast.core.StateNode#copy()
     */
    @Override
    public StateNode copy() {
        Reconciliation copy = new Reconciliation();
        copy.embeddedTree = embeddedTree;
        copy.hostTree = hostTree;
        System.arraycopy(map, 0, copy.map, 0, map.length);
        System.arraycopy(storedMap, 0, copy.storedMap, 0, storedMap.length);
        return copy;
    }

    /* (non-Javadoc)
     * @see beast.core.StateNode#assignTo(beast.core.StateNode)
     */
    @Override
    public void assignTo(StateNode other) {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see beast.core.StateNode#assignFrom(beast.core.StateNode)
     */
    @Override
    public void assignFrom(StateNode other) {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see beast.core.StateNode#assignFromFragile(beast.core.StateNode)
     */
    @Override
    public void assignFromFragile(StateNode other) {
        // TODO Auto-generated method stub

    }

    public String toString() {
        StringBuilder sb = new StringBuilder("");
        for (int i = 0; i < map.length; ++i) {
            if (i > 0) sb.append(" ");
            sb.append(map[i]);
        }
        return sb.toString();
    }
    
    /* (non-Javadoc)
     * @see beast.core.StateNode#fromXML(org.w3c.dom.Node)
     */
    @Override
    public void fromXML(org.w3c.dom.Node node) {
        
        // TODO finish implementation
        String[] tokens = node.getTextContent().split(" ");
        for (int i = 0; i < tokens.length; ++i)
            map[i] = Integer.parseInt(tokens[i]);

    }

    /* (non-Javadoc)
     * @see beast.core.StateNode#scale(double)
     */
    @Override
    public int scale(double fScale) throws Exception {
        // TODO Auto-generated method stub
        return 0;
    }

    /* (non-Javadoc)
     * @see beast.core.StateNode#store()
     */
    @Override
    protected void store() {
        System.arraycopy(map, 0, storedMap, 0, map.length);
    }

    /* (non-Javadoc)
     * @see beast.core.StateNode#restore()
     */
    @Override
    public void restore() {
        System.arraycopy(storedMap, 0, map, 0, map.length);
    }

    public Node getHost(Node embeddedNode) {
        
        if (embeddedNode.getTree() != embeddedTree)
            throw new IllegalArgumentException("embeddedNode must belong to "
                    + "embeddedTree");

        return hostTree.getNode(map[embeddedNode.getNr()]);
        
    }
    
    public void setHost(Node embeddedNode, Node hostNode) {
        
        if (!hasStartedEditing)
            throw new RuntimeException("Have not started editing!");
        
        if (embeddedNode.getTree() != embeddedTree)
            throw new IllegalArgumentException("embeddedNode must belong to "
                    + "embeddedTree");
        
        if (hostNode.getTree() != hostTree)
            throw new IllegalArgumentException("hostNode must belong to "
                    + "hostTree");

        map[embeddedNode.getNr()] = hostNode.getNr();
        
    }
    
    public boolean notarizeTrees(Tree embeddedTree, Tree hostTree) {
        
        return embeddedTree == this.embeddedTree &&
                hostTree == this.hostTree;
    
    }
    
}
