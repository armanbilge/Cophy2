/**
 * ThresholdedCophylogenyModel.java
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

import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.commons.math.util.MathUtils;

import cophy.Utils;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class ThresholdedCophylogenyModel extends EmbeddedTreeDistribution {
    
    public Input<Integer> thresholdInput = new Input<Integer>("threshold",
            "A threshold on the total number of symbionts that may spawn",
            Validate.REQUIRED);
    
    protected int threshold;
        
    public void initAndValidate() throws Exception {
        
        super.initAndValidate();
        threshold = thresholdInput.get();
        
    }

    protected double calculateDensity() {
        
        Tree embeddedTree = embeddedTreeInput.get();
        Tree hostTree = hostTreeInput.get();
        
        SortedMap<Double,Node> hostTreeIterator = new TreeMap<Double,Node>();
        for (Node node : hostTree.getInternalNodes())
            hostTreeIterator.put(node.getHeight(), node);
        return calculateDensity();
        
    }
    
    protected double calculateDensity(double origin, Node embedded) {
        
        return 0.0;
    }
    
    protected DoubleMatrix2D constructRateMatrix(int hostCount) {
        
        double lambda = duplicationRateParameterInput.get().getValue();
        double tau = hostSwitchRateParameterInput.get().getValue();
        double mu = lossRateParameterInput.get().getValue();
        
        int stateCount = MathUtils.pow(threshold + 1, hostCount);
        int[][] decomposedStates = new int[stateCount][];
        for (int i = 0; i < stateCount; ++i)
            decomposedStates[i] = decomposeState(i, hostCount);
        
        DoubleMatrix2D matrix =
                DoubleFactory2D.sparse.make(stateCount, stateCount);
        
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                
                double rate = 0.0;
                int[] state = decomposedStates[i];
                
                if (i == j) {
                    if (i < stateCount - 1) rate += lambda + tau;
                    rate += mu;
                    rate *= -Utils.sum(state);
                } else {
                    
                    int eventLocus = locateBirthEvent(state,
                            decomposedStates[j]);
                    if (eventLocus != -1) {
                        rate += state[eventLocus] * lambda;
                        rate += (Utils.sum(state) - state[eventLocus]) /
                                (hostCount - 1) * tau;
                    }
                    
                    eventLocus = locateBirthEvent(decomposedStates[j], state);
                    if (eventLocus != -1)
                        rate += state[eventLocus] * mu;
                    
                }
                
                if (rate != 0.0) matrix.set(i, j, rate);
                
            }
        }
        
        return matrix;
    }
    
    protected final int[] decomposeState(int state, int hostCount) {
        
        int base = threshold + 1;
        int[] decomposedState = new int[hostCount];
        for (int i = 0; i < hostCount; ++i) {
            decomposedState[i] = state % base;
            state /= base;
        }
        return decomposedState;
        
    }
    
    
    protected final int locateBirthEvent(int[] s0, int[] s1) {
        
        int eventLocus = -1;
        
        for (int i = 0; i < s0.length; ++i) {
                        
            int change = s1[i] - s0[i];
            
            if (change > 1 || change < 0) {
                return -1;
            } else if (change == 1) {
                if (eventLocus == - 1)
                    eventLocus = i;
                else
                    return -1;
            }
            
        }
                
        return eventLocus;
        
    }

    @Override
    protected double getOverallRate(Node embedded, Node host) {
        
        BranchRateModel branchRateModel = branchRateModelInput.get();
        return branchRateModel.getRateForBranch(embedded);
        
    }
    
}
