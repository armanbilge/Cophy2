/**
 * HostSwitchOperator.java
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

import java.util.List;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import cophy.Reconciliation;
import cophy.Utils;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class HostSwitchOperator extends Operator {

    public Input<Tree> embeddedTreeInput = new Input<Tree>("embeddedTree",
            "The tree embedded within the host tree.",
            Input.Validate.REQUIRED);
    public Input<Tree> hostTreeInput = new Input<Tree>("hostTree",
            "The host, or container, tree.", Input.Validate.REQUIRED);
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
        Tree hostTree = hostTreeInput.get(this);
        Reconciliation reconciliation = reconciliationInput.get(this);
        double originHeight = originHeightParameterInput.get(this).getValue();
        
        Node embeddedNode = embeddedTree.getInternalNodes()
                .get(Randomizer.nextInt(embeddedTree.getInternalNodeCount()));
        
        if (embeddedNode.getHeight()
                == reconciliation.getHost(embeddedNode).getHeight())
             return Double.NEGATIVE_INFINITY;
        
        double lower = Math.max(embeddedNode.getLeft().getHeight(),
                embeddedNode.getRight().getHeight());
        double upper = embeddedNode.isRoot() ? originHeight :
                    embeddedNode.getParent().getHeight();
        double range = upper - lower;
        
        double newHeight = Randomizer.nextDouble() * range + lower;
        
        List<Node> potentialHosts =
                Utils.getLineagesAtHeight(hostTree, newHeight);
        Node newHost =
                potentialHosts.get(Randomizer.nextInt(potentialHosts.size()));
        
        embeddedNode.setHeight(newHeight);
        reconciliation.setHost(embeddedNode, newHost);
        
        int inversePotentialHosts = Utils.getLineageCountAtHeight(hostTree,
                embeddedNode.getHeight(), false);
        
        return Math.log(potentialHosts.size())
                - Math.log(inversePotentialHosts);
    }

}
