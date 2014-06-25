/**
 * CospeciationOperator.java
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

package cophy.operator;

import cophy.Reconciliation;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.MachineAccuracy;
import beast.util.Randomizer;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class CospeciationOperator extends Operator {

    public Input<Tree> embeddedTreeInput = new Input<Tree>("embeddedTree",
            "The tree embedded within the host tree.",
            Input.Validate.REQUIRED);
    public Input<Reconciliation> reconciliationInput =
            new Input<Reconciliation>("reconciliation", "A mapping of nodes of"
                    + " the embedded tree to nodes of its host tree",
                    Input.Validate.REQUIRED);
    public Input<RealParameter> originHeightParameterInput =
            new Input<RealParameter>("originHeight", "Parameter"
                    + " specifying the origin of the embedded process",
                    Input.Validate.REQUIRED);
    
    public void initAndValidate() {}
    
    @Override
    public double proposal() {
        
        Tree embeddedTree = embeddedTreeInput.get(this);
        Reconciliation reconciliation = reconciliationInput.get(this);
        double originHeight = originHeightParameterInput.get(this).getValue();
        
        Node embeddedNode = embeddedTree.getInternalNodes()
                .get(Randomizer.nextInt(embeddedTree.getInternalNodeCount()));
        Node hostNode = reconciliation.getHost(embeddedNode);
        if (hostNode.isLeaf()) return Double.NEGATIVE_INFINITY;
        double hostHeight = hostNode.getHeight();
        
        double child0Height = embeddedNode.getLeft().getHeight();
        double child1Height = embeddedNode.getRight().getHeight();
        if (!embeddedNode.isLeaf()
                && hostHeight <= Math.max(child0Height, child1Height))
            return Double.NEGATIVE_INFINITY;
        
        double upperHeightEmbedded = embeddedNode.isRoot() ? originHeight :
            embeddedNode.getParent().getHeight();
        double upperHeightHost = hostNode.isRoot() ? Double.POSITIVE_INFINITY :
            hostNode.getParent().getHeight();
        double upperHeight = Math.min(upperHeightEmbedded, upperHeightHost);
        double range = upperHeight - hostHeight;
        
        double newHeight, hastingsRatio;
        if (MachineAccuracy.same(embeddedNode.getHeight(), hostHeight)) {
            newHeight = Randomizer.nextDouble() * range + hostHeight;
            hastingsRatio = range;
        }
        else {
            newHeight = hostHeight;
            hastingsRatio = 1 / range;
        }
        
        embeddedNode.setHeight(newHeight);
        return Math.log(hastingsRatio);
        
    }

}
