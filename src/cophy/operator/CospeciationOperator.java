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

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import cophy.Reconciliation;

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
        
        double leftChildHeight = embeddedNode.getLeft().getHeight();
        double rightChildHeight = embeddedNode.getRight().getHeight();
        if (hostHeight <= Math.max(leftChildHeight, rightChildHeight))
            return Double.NEGATIVE_INFINITY;
        
        double upperHeightEmbedded = embeddedNode.isRoot() ? originHeight :
            embeddedNode.getParent().getHeight();
        double upperHeightHost = hostNode.isRoot() ? Double.POSITIVE_INFINITY :
            hostNode.getParent().getHeight();
        double upperHeight = Math.min(upperHeightEmbedded, upperHeightHost);
        double range = upperHeight - hostHeight;
        
        double newHeight;
        double logHastingsRatio = Math.log(range);
        if (embeddedNode.getHeight() == hostHeight) {
            newHeight = Randomizer.nextDouble() * range + hostHeight;
        } else {
            newHeight = hostHeight;
            logHastingsRatio *= -1;
        }
        
        embeddedNode.setHeight(newHeight);
        return logHastingsRatio;
    }

}
