/**
 * LeafHostSwitchOperator.java
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

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import cophy.Reconciliation;


/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class LeafHostSwitchOperator extends HostSwitchOperator {

    @Override
    public double proposal() {

        Tree embeddedTree = embeddedTreeInput.get(this);
        Tree hostTree = hostTreeInput.get(this);
        Reconciliation reconciliation = reconciliationInput.get(this);
        
        Node embeddedNode = embeddedTree.getExternalNodes()
                .get(Randomizer.nextInt(embeddedTree.getLeafNodeCount()));
        Node hostNode = reconciliation.getHost(embeddedNode);
        Node newHostNode = hostTree.getExternalNodes()
                .get(Randomizer.nextInt(hostTree.getLeafNodeCount()));
        
        if (hostNode.equals(newHostNode)) return Double.NEGATIVE_INFINITY;
        
        reconciliation.setHost(embeddedNode, newHostNode);
        
        return 0.0;
    }

}
