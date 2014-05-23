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

import cophy.Reconciliation;
import cophy.Utils;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.Parameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class HostSwitchOperator extends Operator {

    public Input<Tree> embeddedTreeInput = new Input<Tree>("embedded tree",
            "The tree embedded within the host tree.",
            Input.Validate.REQUIRED);
    public Input<Tree> hostTreeInput = new Input<Tree>("host tree",
            "The host, or container, tree.", Input.Validate.REQUIRED);
    public Input<Reconciliation> reconciliationInput =
            new Input<Reconciliation>("reconciliation", "A mapping of nodes of"
                    + " the embedded tree to nodes of its host tree",
                    Input.Validate.REQUIRED);
    public Input<Parameter<Double>> originHeightParameterInput =
            new Input<Parameter<Double>>("origin height", "Parameter"
                    + " specifying the origin of the embedded process",
                    Input.Validate.REQUIRED);
    
    protected final Tree embeddedTree;
    protected final Tree hostTree;
    protected final Reconciliation reconciliation;
    protected final Parameter<Double> originHeightParameter;
    
    /**
     * 
     */
    public HostSwitchOperator(Tree embeddedTree, Tree hostTree,
            Reconciliation reconciliation, Parameter<Double> originHeight) {
        
        this.embeddedTree = embeddedTree;
        this.hostTree = hostTree;
        this.reconciliation = reconciliation;
        this.originHeightParameter = originHeight;
        
    }

    @Override
    public double proposal() {
        
        Node embeddedNode = embeddedTree.getInternalNodes()
                .get(Randomizer.nextInt(embeddedTree.getInternalNodeCount()));
        double lower = Math.max(embeddedNode.getChild(0).getHeight(),
                embeddedNode.getChild(1).getHeight());
        double upper = embeddedNode.isRoot() ?
                originHeightParameter.getArrayValue() :
                    embeddedNode.getParent().getHeight();
        List<Node> potentialNewHosts =
                Utils.getLineagesInHeightRange(hostTree, lower, upper);
        Node proposedHostNode = potentialNewHosts.get(
                Randomizer.nextInt(potentialNewHosts.size()));
        Node currentHostNode = reconciliation.getHost(embeddedNode);
        if (proposedHostNode == currentHostNode) {
            reject(-2);
            return 0.0;
        }
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
