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
    
    protected double calculateDensity(double originHeight, int startState,
            Node embedded, TreeMap<Double,Node> hostSpeciations,
            DoubleMatrix2D[] matrices) {
        
        double L = 1.0;
        
        Tree hostTree = hostTreeInput.get();
        Reconciliation reconciliation = reconciliationInput.get();
        
        BranchRateModel branchRateModel = branchRateModelInput.get();
        double rate = branchRateModel.getRateForBranch(embedded);
        
        double embeddedHeight = embedded.getHeight();
        
        Node host = reconciliation.getHost(embedded);
        if (!Utils.lineageExistedAtHeight(host, embeddedHeight))
            return Double.NEGATIVE_INFINITY;
        
        double startHeight = originHeight;
        
        Map<Node,Integer> hostNodes2Bins = new HashMap<Node,Integer>();
        List<Node> startHostLineages = Utils.getLineagesAtHeight(hostTree,
                startHeight, true);
        int hostCount = startHostLineages.size();
        for (int i = 0; i < hostCount; ++i)
            hostNodes2Bins.put(startHostLineages.get(i), i);
        
        NavigableSet<Double> hostSpeciationSet = hostSpeciations.subMap(
                embeddedHeight, true, startHeight, false).descendingKeySet();
        
        DoubleMatrix1D startDensity =
                DoubleFactory1D.dense.make(getStateCount(hostCount));
        startDensity.setQuick(startState, 1.0);
        
        int speciatedBin;
        DoubleMatrix1D endDensity = DoubleFactory1D.sparse.make(hostCount);
        
        for (double speciationHeight : hostSpeciationSet) {
            
            DoubleMatrix2D matrix = matrices[hostCount-1];
            double t = (startHeight - speciationHeight)
                    * rate;
            endDensity = AMH11.expmv(t, matrix, startDensity);
            startHeight = speciationHeight;
            speciatedBin = hostNodes2Bins.get(
                    hostSpeciations.get(speciationHeight));
            hostNodes2Bins.clear();
            ++hostCount;
            startHostLineages =
                    Utils.getLineagesAtHeight(hostTree, startHeight, true);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            int[] map = mapNewStatesToOld(hostCount, speciatedBin);
            startDensity = DoubleFactory1D.dense.make(getStateCount(hostCount));
            for (int i = 0; i < map.length; ++i) {
                double density = map[i] != -1 ? endDensity.getQuick(map[i]) : 0;
                startDensity.setQuick(i, density);
            }
            
        }
        
        Node embeddedLeft = embedded.getLeft();
        Node embeddedRight = embedded.getRight();
        int[] state = new int[hostCount];
        
        if (MachineAccuracy.same(embeddedHeight, 0.0)) {
            
            int hostBin = hostNodes2Bins.get(host);
            state[hostBin] = 1;
            return endDensity.getQuick(compressState(state));
            
        } else if (MachineAccuracy.same(embeddedHeight, host.getHeight())) {

            int hostBin = hostNodes2Bins.get(host.getLeft());
            state[hostBin] = 1;
            double p1 = calculateDensity(embeddedHeight, compressState(state),
                    embeddedLeft, hostSpeciations, matrices);
            double p2 = calculateDensity(embeddedHeight, compressState(state),
                    embeddedRight, hostSpeciations, matrices);
            state[hostBin] = 0;

            hostBin = hostNodes2Bins.get(host.getRight());
            state[hostBin] = 1;
            p1 *= calculateDensity(embeddedHeight, compressState(state),
                    embeddedRight, hostSpeciations, matrices);
            p2 *= calculateDensity(embeddedHeight, compressState(state),
                    embeddedLeft, hostSpeciations, matrices); 
            state[hostBin] = 0;
            
            L *= p1 + p2;
            
            int[] map = mapStatesToOneLess(hostCount-1, hostBin);
            int[] map1 = mapNewStatesToOld(hostCount, hostBin);
            
            startDensity = DoubleFactory1D.dense.make(hostCount);
            for (int i = 0; i < hostCount; ++i) {
                double density = map1[i] != -1 && map[map1[i]] != -1 ?
                        endDensity.getQuick(map[map1[i]]) : 0;
                startDensity.setQuick(i, density);
            }
            
            L *= calculateDensity(embeddedHeight, startDensity, hostSpeciations,
                    rate, matrices);
            
            
        } else {
            
            --hostCount;
            double lambda = duplicationRateParameterInput.get().getValue();
            double tau = hostSwitchRateParameterInput.get().getValue();
            double pDuplication = lambda / (lambda + tau);
            L *= bdSpeciationDensity(embeddedHeight, originHeight);
            int hostBin = hostNodes2Bins.get(host);
                        
            Node integrated = embeddedLeft;
            state[hostBin] = 1;
            double sum1 = pDuplication * calculateDensity(embeddedHeight,
                    compressState(state), integrated, hostSpeciations,
                    matrices);
            state[hostBin] = 0;
            double sum2 = 0.0;
            for (int i = 0; i < hostCount; ++i) {
                if (i != hostBin) {
                    state[i] = 1;
                    sum2 += calculateDensity(embeddedHeight,
                            compressState(state), integrated, hostSpeciations,
                            matrices);
                    state[i] = 0;
                }
            }
            sum1 += (1 - pDuplication) * sum2 / (hostCount - 1);
            
            state[hostBin] = 1;
            double sum = sum1 * calculateDensity(embeddedHeight,
                    compressState(state), embeddedRight, hostSpeciations,
                    matrices);
            state[hostBin] = 0;

            integrated = embeddedRight;
            state[hostBin] = 1;
            sum1 = pDuplication * calculateDensity(embeddedHeight,
                    compressState(state), integrated, hostSpeciations,
                    matrices);
            state[hostBin] = 0;
            sum2 = 0.0;
            for (int i = 0; i < hostCount; ++i) {
                if (i != hostBin) {
                    state[i] = 1;
                    sum2 += calculateDensity(embeddedHeight,
                            compressState(state), integrated, hostSpeciations,
                            matrices);
                    state[i] = 0;
                }
            }
            sum1 += (1 - pDuplication) * sum2 / (hostCount - 1);

            state[hostBin] = 1;
            sum += sum1 * calculateDensity(embeddedHeight, compressState(state),
                    embeddedLeft, hostSpeciations, matrices);
            state[hostBin] = 0;
                        
            L *= sum;
            
            int[] map = mapStatesToOneLess(hostCount, hostBin);
            startDensity = DoubleFactory1D.dense.make(hostCount);
            for (int i = 0; i < hostCount; ++i) {
                double density = map[i] != -1 ? endDensity.getQuick(map[i]) : 0;
                startDensity.setQuick(i, density);
            }
            
            L *= calculateDensity(embeddedHeight, startDensity, hostSpeciations,
                    rate, matrices);
            
        }
                
        return L;
        
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
        DoubleMatrix1D endDensity = DoubleFactory1D.sparse.make(hostCount);
        
        for (double speciationHeight : hostSpeciationSet) {
            
            DoubleMatrix2D matrix = matrices[hostCount-1];
            double t = (startHeight - speciationHeight) * rate;
            endDensity = AMH11.expmv(t, matrix, startDensity);
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
        
        return endDensity.get(0);
    }
    
    protected double bdSpeciationDensity(double s, double t) {
        
        double lambda = duplicationRateParameterInput.get().getValue()
            + hostSwitchRateParameterInput.get().getValue();
        double mu = lossRateParameterInput.get().getValue();
        
        double a = lambda - mu;
        double b = Math.exp(-a * s);
        double c = lambda - mu * b;
        double d = Math.exp(-a * t);
        
        return a * a * b / (c * c) * (lambda - mu * d) / (1 - d);
        
    }
    
    protected DoubleMatrix2D constructRateMatrix(int hostCount) {
        
        double lambda = duplicationRateParameterInput.get().getValue();
        double tau = hostSwitchRateParameterInput.get().getValue();
        double mu = lossRateParameterInput.get().getValue();
        
        int stateCount = getStateCount(hostCount);
        int[][] decompressedStates = new int[stateCount][];
        for (int i = 0; i < stateCount; ++i)
            decompressedStates[i] = decompressState(i, hostCount);
        
        DoubleMatrix2D matrix =
                DoubleFactory2D.sparse.make(stateCount, stateCount);
        
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                
                double rate = 0.0;
                int[] state = decompressedStates[i];
                
                if (i == j) {
                    if (i < stateCount - 1) rate += lambda + tau;
                    rate += mu;
                    rate *= -Utils.sum(state);
                } else {
                    
                    int eventLocus = locateBirthEvent(state,
                            decompressedStates[j]);
                    if (eventLocus != -1) {
                        rate += state[eventLocus] * lambda;
                        rate += (Utils.sum(state) - state[eventLocus]) /
                                (hostCount - 1) * tau;
                    }
                    
                    eventLocus = locateBirthEvent(decompressedStates[j], state);
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
    
    protected final int compressState(int[] decompressedState) {
        
        int base = threshold + 1;
        int state = 0;
        for (int i = 0, j = 1; i < decompressedState.length; ++i, j *= base)
            state += decompressedState[i] * j;
        return state;
        
    }
    
    protected final int[] decompressState(int state, int hostCount) {
        int[] decompressedState = new int[hostCount];
        decompressState(state, decompressedState);
        return decompressedState;
    }
    
    protected final void decompressState(int state, int[] decompressedState) {
        
        int base = threshold + 1;
        for (int i = 0; i < decompressedState.length; ++i) {
            decompressedState[i] = state % base;
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
        
        int[] s0Decompressed = decompressState(s0, hostCount);
        int[] s1Decompressed = decompressState(s1, hostCount);
        for (int i = 0; i < hostCount; ++i) {
            if (s0Decompressed[1] < s1Decompressed[0]) return false;
        }
        return true;
        
    }

    protected final int[] mapNewStatesToOld(int hostCount, int speciatedBin) {
        
        int newStateCount = getStateCount(hostCount);
        int[] decompressedState = new int[hostCount];
        int[] map = new int[getStateCount(hostCount)];
        int[] oldState = new int[hostCount-1];
        for (int state = 0; state < newStateCount; ++state) {
            boolean valid = true;
            decompressState(state, decompressedState);
            for (int i = 0; i < hostCount; ++i) {
                int j = i > speciatedBin ? i : i-1;
                oldState[j] = decompressedState[i];
                if (i == speciatedBin) {
                    if (decompressedState[i] != decompressedState[i+1]) {
                        valid = false;
                        break;
                    }
                    ++i;
                }
            }
            if (valid)
                map[state] = compressState(oldState);
            else
                map[state] = -1;
        }
        return map;
        
    }
    
    protected final int[] mapStatesToOneLess(int hostCount, int speciatedBin) {
        
        int newStateCount = getStateCount(hostCount);
        int[] decompressedState = new int[hostCount];
        int[] map = new int[getStateCount(hostCount)];
        for (int state = 0; state < newStateCount; ++state) {
            decompressState(state, decompressedState);
            if (decompressedState[speciatedBin] < 1) {
                map[state] = -1;
            } else {
                decompressedState[speciatedBin] -= 1;
                map[state] = compressState(decompressedState);
            }
        }
        return map;
        
    }
    
    @Override
    protected double getOverallRate(Node embedded, Node host) {
        
        BranchRateModel branchRateModel = branchRateModelInput.get();
        return branchRateModel.getRateForBranch(embedded);
        
    }
    
}
