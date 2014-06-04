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
        double lower = Math.max(embeddedNode.getLeft().getHeight(),
                embeddedNode.getRight().getHeight());
        double upper = embeddedNode.isRoot() ? originHeight :
                    embeddedNode.getParent().getHeight();
        if (upper < lower) return Double.NEGATIVE_INFINITY;
        List<Node> potentialNewHosts =
                Utils.getLineagesInHeightRange(hostTree, lower, upper);
        Node proposedHostNode = potentialNewHosts.get(
                Randomizer.nextInt(potentialNewHosts.size()));
        Node currentHostNode = reconciliation.getHost(embeddedNode);
        if (proposedHostNode.equals(currentHostNode))
            return Double.NEGATIVE_INFINITY;
        
        double hastingsRatio = 1.0;
        if (!Utils.lineageExistedAtHeight(proposedHostNode,
                embeddedNode.getHeight())) {
            double hastingsLower =
                    Math.max(lower, proposedHostNode.getHeight());
            double hastingsUpper = Math.min(upper,
                    proposedHostNode.isRoot() ? Double.POSITIVE_INFINITY :
                        proposedHostNode.getParent().getHeight());
            embeddedTree.startEditing(this);
            embeddedNode.setHeight(Randomizer.nextDouble()
                    * (upper - lower) + lower);
            double inverseHastingsLower =
                    Math.max(lower, currentHostNode.getHeight());
            double inverseHastingsUpper = Math.min(upper,
                    currentHostNode.isRoot() ? Double.POSITIVE_INFINITY :
                        currentHostNode.getParent().getHeight());
            hastingsRatio = (inverseHastingsUpper - inverseHastingsLower) /
                    (hastingsUpper - hastingsLower);
        }
        return Math.log(hastingsRatio);
    }

}
