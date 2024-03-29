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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class Reconciliation extends IntegerParameter implements Function {

    public Input<Tree> embeddedTreeInput = new Input<Tree>("embeddedTree", 
            "The embedded, or contained, tree.", Validate.REQUIRED);
    public Input<Tree> hostTreeInput = new Input<Tree>("hostTree",
            "The tree hosting the embedded tree.", Validate.REQUIRED);
    public Input<TraitSet> traitSetInput = new Input<TraitSet>("hostTrait",
            "The set of host traits for the embedded taxa.", Validate.REQUIRED);
    public Input<RealParameter> originHeightParameterInput =
            new Input<RealParameter>("originHeight", "The height of origin" +
                    " for the embedded tree", Validate.REQUIRED);
        
    @Override
    public void initAndValidate() throws Exception {
        
        super.initAndValidate();
        Tree embeddedTree = embeddedTreeInput.get();
        setDimension(embeddedTree.getNodeCount());
        Tree hostTree = hostTreeInput.get();
        setBounds(0, hostTree.getNodeCount() - 1);
        TraitSet hostTraitSet = traitSetInput.get();
        
        RealParameter originHeightParameter = originHeightParameterInput.get();
        
        Map<String,Node> hostTaxon2Node = new HashMap<String,Node>();
        for (Node node : hostTree.getExternalNodes()) {
            hostTaxon2Node.put(hostTree.getTaxonId(node), node);
        }        
        
        for (Node embeddedNode : embeddedTree.getExternalNodes()) {
            String hostTaxon =
                    hostTraitSet.getStringValue(embeddedNode.getNr());
            Node hostNode = hostTaxon2Node.get(hostTaxon);
            setHostInit(embeddedNode, hostNode);
        }
        
        for (Node embeddedNode : embeddedTree.getInternalNodes()) {
            List<Node> potentialHosts = Util.getLineagesAtHeight(hostTree,
                    embeddedNode.getHeight());
            int r = Randomizer.nextInt(potentialHosts.size());
            assert(Util.lineageExistedAtHeight(potentialHosts.get(r), embeddedNode.getHeight()));
            setHostInit(embeddedNode, potentialHosts.get(r));
        }
        
        originHeightParameter.setBounds(0.0, Double.POSITIVE_INFINITY);
        
    }
    
    public Node getHost(Node embeddedNode) {
        int hostNodeNr = getNativeValue(embeddedNode.getNr());
        return hostTreeInput.get().getNode(hostNodeNr);
    }
    
    public void setHost(Node embeddedNode, Node hostNode) {
        setValue(embeddedNode.getNr(), hostNode.getNr());
    }

    private void setHostInit(Node embeddedNode, Node hostNode) {
        
        int iParam = embeddedNode.getNr();
        int fValue = hostNode.getNr();
        values[iParam] = fValue;
        m_bIsDirty[iParam] = true;
        m_nLastDirty = iParam;
        
    }
        
}
