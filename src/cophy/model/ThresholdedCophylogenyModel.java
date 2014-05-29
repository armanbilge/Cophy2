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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeMap;

import org.apache.commons.math.util.MathUtils;

import amh11.AMH11;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.MachineAccuracy;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cophy.Reconciliation;
import cophy.Utils;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class ThresholdedCophylogenyModel extends EmbeddedTreeDistribution {
    
    public Input<Integer> thresholdInput = new Input<Integer>("threshold",
            "A threshold on the total number of symbionts that may spawn",
            Validate.REQUIRED);
    
    protected int threshold;
        
    public void initAndValidate() {
        
        super.initAndValidate();
        threshold = thresholdInput.get();
        
    }

    protected double calculateDensity() {
        
        Tree embeddedTree = embeddedTreeInput.get();
        Tree hostTree = hostTreeInput.get();
        double originHeight = originHeightParameterInput.get().getValue();
        
        if (originHeight < embeddedTree.getRoot().getHeight()
                || originHeight < hostTree.getRoot().getHeight())
            return Double.NEGATIVE_INFINITY;
        
        TreeMap<Double,Node> hostSpeciations = new TreeMap<Double,Node>();
        for (Node node : hostTree.getInternalNodes())
            hostSpeciations.put(node.getHeight(), node);
        
        DoubleMatrix2D[] matrices =
                new DoubleMatrix2D[hostTree.getLeafNodeCount()-1];
        
        for (int i = 0; i < matrices.length; ++i) {
            matrices[i] = constructRateMatrix(i+1);
        }
        
        return calculateDensity(originHeight, 1, embeddedTree.getRoot(),
                hostSpeciations, matrices);
        
    }
    
    protected double calculateDensity(double startHeight, int startState,
            Node embedded, TreeMap<Double,Node> hostSpeciations,
            DoubleMatrix2D[] matrices) {
        
        Tree hostTree = hostTreeInput.get();
        Reconciliation reconciliation = reconciliationInput.get();
        
        BranchRateModel branchRateModel = branchRateModelInput.get();
        
        Node host = reconciliation.getHost(embedded);
        if (!Utils.lineageExistedAtHeight(host, embedded.getHeight()))
            return Double.NEGATIVE_INFINITY;
        
        Map<Node,Integer> hostNodes2Bins = new HashMap<Node,Integer>();
        List<Node> startHostLineages = Utils.getLineagesAtHeight(hostTree,
                startHeight);
        int hostCount = startHostLineages.size();
        for (int i = 0; i < hostCount; ++i)
            hostNodes2Bins.put(startHostLineages.get(i), i);
        
        NavigableSet<Double> hostSpeciationSet = hostSpeciations.subMap(
                embedded.getHeight(), true, startHeight, false)
                .descendingKeySet();
        
        DoubleMatrix1D startDensity =
                DoubleFactory1D.dense.make(getStateCount(hostCount));
        startDensity.setQuick(startState, 1.0);
        
        int speciatedBin;
        
        for (double speciationHeight : hostSpeciationSet) {
            
            DoubleMatrix2D matrix = matrices[hostCount-1];
            double t = (startHeight - speciationHeight)
                    * branchRateModel.getRateForBranch(embedded);
            DoubleMatrix1D endDensity = AMH11.expmv(t, matrix, startDensity);
            startHeight = speciationHeight;
            speciatedBin = hostNodes2Bins.get(
                    hostSpeciations.get(speciationHeight));
            hostNodes2Bins.clear();
            ++hostCount;
            startHostLineages =
                    Utils.getLineagesAtHeight(hostTree, startHeight);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            int[] map = mapNewStatesToOld(hostCount, speciatedBin);
            startDensity = DoubleFactory1D.dense.make(getStateCount(hostCount));
            for (int i = 0; i < map.length; ++i) {
                double density = map[i] != -1 ? endDensity.getQuick(map[i]) : 0;
                startDensity.setQuick(i, density);
            }
            
        }
        
        if (MachineAccuracy.same(embedded.getHeight(), host.getHeight())) {
            // TODO
        } else {
            // TODO
        }
        
        return 0.0;
        
    }
    
    protected double calculateDensity(double startHeight,
            DoubleMatrix1D startDensity, TreeMap<Double,Node> hostSpeciations,
            double rate, DoubleMatrix2D[] matrices) {
        
        Tree hostTree = hostTreeInput.get();
                        
        Map<Node,Integer> hostNodes2Bins = new HashMap<Node,Integer>();
        List<Node> startHostLineages = Utils.getLineagesAtHeight(hostTree,
                startHeight);
        int hostCount = startHostLineages.size();
        for (int i = 0; i < hostCount; ++i)
            hostNodes2Bins.put(startHostLineages.get(i), i);
        
        NavigableSet<Double> hostSpeciationSet = hostSpeciations.subMap(0.0,
                true, startHeight, false).descendingKeySet();
                
        int speciatedBin;
        
        for (double speciationHeight : hostSpeciationSet) {
            
            DoubleMatrix2D matrix = matrices[hostCount-1];
            double t = (startHeight - speciationHeight) * rate;
            DoubleMatrix1D endDensity = AMH11.expmv(t, matrix, startDensity);
            startHeight = speciationHeight;
            speciatedBin = hostNodes2Bins.get(
                    hostSpeciations.get(speciationHeight));
            hostNodes2Bins.clear();
            ++hostCount;
            startHostLineages =
                    Utils.getLineagesAtHeight(hostTree, startHeight);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            int[] map = mapNewStatesToOld(hostCount, speciatedBin);
            startDensity = DoubleFactory1D.dense.make(getStateCount(hostCount));
            for (int i = 0; i < map.length; ++i) {
                double density = map[i] != -1 ? endDensity.getQuick(map[i]) : 0;
                startDensity.setQuick(i, density);
            }
            
        }
        
        return Math.log(startDensity.get(0));
    }
    
    protected DoubleMatrix2D constructRateMatrix(int hostCount) {
        
        double lambda = duplicationRateParameterInput.get().getValue();
        double tau = hostSwitchRateParameterInput.get().getValue();
        double mu = lossRateParameterInput.get().getValue();
        
        int stateCount = getStateCount(hostCount);
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
    
    protected final int getStateCount(int hostCount) {
        return MathUtils.pow(threshold + 1, hostCount);
    }
    
    protected final int composeState(int[] decomposedState) {
        
        int base = threshold + 1;
        int state = 0;
        for (int i = 0, j = 1; i < decomposedState.length; ++i, j *= base)
            state += decomposedState[i] * j;
        return state;
        
    }
    
    protected final int[] decomposeState(int state, int hostCount) {
        int[] decomposedState = new int[hostCount];
        decomposeState(state, decomposedState);
        return decomposedState;
    }
    
    protected final void decomposeState(int state, int[] decomposedState) {
        
        int base = threshold + 1;
        for (int i = 0; i < decomposedState.length; ++i) {
            decomposedState[i] = state % base;
            state /= base;
        }
        
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
    
    protected final boolean statesAreCompatible(int s0, int s1, int hostCount) {
        
        int[] s0Decomposed = decomposeState(s0, hostCount);
        int[] s1Decomposed = decomposeState(s1, hostCount);
        for (int i = 0; i < hostCount; ++i) {
            if (s0Decomposed[1] < s1Decomposed[0]) return false;
        }
        return true;
    }

    protected final int[] mapNewStatesToOld(int hostCount, int speciatedBin) {
        
        int newStateCount = getStateCount(hostCount);
        int[] decomposedState = new int[hostCount];
        int[] map = new int[getStateCount(hostCount)];
        int[] oldState = new int[hostCount-1];
        for (int state = 0; state < newStateCount; ++state) {
            boolean valid = true;
            decomposeState(state, decomposedState);
            for (int i = 0; i < hostCount; ++i) {
                int j = i > speciatedBin ? i : i-1;
                oldState[j] = decomposedState[i];
                if (i == speciatedBin) {
                    if (decomposedState[i] != decomposedState[i+1]) {
                        valid = false;
                        break;
                    }
                    ++i;
                }
            }
            if (valid)
                map[state] = composeState(oldState);
            else
                map[state] = -1;
        }
        
        return null;
        
    }
    
    @Override
    protected double getOverallRate(Node embedded, Node host) {
        
        BranchRateModel branchRateModel = branchRateModelInput.get();
        return branchRateModel.getRateForBranch(embedded);
        
    }
    
}
