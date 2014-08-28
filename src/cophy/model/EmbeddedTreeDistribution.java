/**
 * EmbeddedTreeDistribution.java
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

package cophy.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.CalculationNode;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import cophy.Reconciliation;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public abstract class EmbeddedTreeDistribution extends Distribution {

    public Input<Tree> embeddedTreeInput = new Input<Tree>("embeddedTree", 
            "The embedded, or contained, tree.", Validate.REQUIRED);
    public Input<Reconciliation> reconciliationInput =
            new Input<Reconciliation>("reconciliation",
                    "The tree-tree mapping.", Validate.REQUIRED);
    public Input<Tree> hostTreeInput = new Input<Tree>("hostTree",
            "The tree hosting the embedded tree.", Validate.REQUIRED);
    public Input<RealParameter> rateParameterInput =
            new Input<RealParameter>("rate", "The overall event rate.",
                    Validate.REQUIRED);
    public Input<RealParameter> originHeightParameterInput =
            new Input<RealParameter>("originHeight", "The height of origin" +
                    " for the embedded tree", Validate.REQUIRED);
    public Input<RealParameter> duplicationRateParameterInput =
            new Input<RealParameter>("duplicationRate",
                    "The model duplication rate", Validate.REQUIRED);
    public Input<RealParameter> hostSwitchRateParameterInput = 
            new Input<RealParameter>("hostSwitchRate",
                    "The model host-switch rate", Validate.REQUIRED);
    public Input<RealParameter> lossRateParameterInput = 
            new Input<RealParameter>("lossRate", "The model loss rate",
                    Validate.REQUIRED);
    
    public void initAndValidate() {
        duplicationRateParameterInput.get()
            .setBounds(0.0, Double.POSITIVE_INFINITY);
        hostSwitchRateParameterInput.get()
            .setBounds(0.0, Double.POSITIVE_INFINITY);
        lossRateParameterInput.get()
            .setBounds(0.0, Double.POSITIVE_INFINITY);
    }
    
    @Override
    public List<String> getArguments() {
        
        List<String> arguments = new ArrayList<String>();
        arguments.add(embeddedTreeInput.get().getID());
        arguments.add(reconciliationInput.get().getID());
        return arguments;
        
    }

    @Override
    public List<String> getConditions() {
        
        List<String> conditions = new ArrayList<String>();
        conditions.add(hostTreeInput.get().getID());
        conditions.add(((CalculationNode) rateParameterInput.get()).getID());
        conditions.add(originHeightParameterInput.get().getID());
        conditions.add(duplicationRateParameterInput.get().getID());
        conditions.add(hostSwitchRateParameterInput.get().getID());
        conditions.add(lossRateParameterInput.get().getID());
        return conditions;
        
    }

    @Override
    public double calculateLogP() throws Exception {
        if (isDirtyCalculation()) logP = calculateLogDensity();
        return logP;
    }
    
    protected abstract double calculateLogDensity();
    
    @Override
    public void sample(State state, Random random) {
        throw new RuntimeException("Not implemented yet!");
    }

}
