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

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class Reconciliation extends CalculationNode {

    public Input<Tree> embeddedTreeInput = new Input<Tree>("embeddedTree", 
            "The embedded, or contained, tree.", Validate.REQUIRED);
    public Input<Tree> hostTreeInput = new Input<Tree>("hostTree",
            "The tree hosting the embedded tree.", Validate.REQUIRED);
    public Input<IntegerParameter> mapInput = new Input<IntegerParameter>(
            "map", "An integer parameter mapping embedded nodes to host nodes",
            Validate.REQUIRED);
        
    @Override
    public void initAndValidate() throws Exception {
        
        super.initAndValidate();
        IntegerParameter map = mapInput.get();
        map.setDimension(embeddedTreeInput.get().getNodeCount());
        map.setLower(0);
        map.setUpper(hostTreeInput.get().getNodeCount() - 1);
        
    }
    
    public Node getHost(Node embeddedNode) {

        int hostNodeNr = mapInput.get().getNativeValue(embeddedNode.getNr());
        return hostTreeInput.get().getNode(hostNodeNr);
        
    }
    
    public void setHost(Node embeddedNode, Node hostNode) {
        
        mapInput.get().setValue(embeddedNode.getNr(), hostNode.getNr());
        
    }
    
}
