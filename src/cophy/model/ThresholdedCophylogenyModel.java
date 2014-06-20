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

    @Override
    protected double calculateLogDensity() {
        
        Tree embeddedTree = embeddedTreeInput.get();
        Tree hostTree = hostTreeInput.get();
        double originHeight = originHeightParameterInput.get().getValue();
        
        if (originHeight < embeddedTree.getRoot().getHeight()
                || originHeight < hostTree.getRoot().getHeight())
            return Double.NEGATIVE_INFINITY;
        
        TreeMap<Double,Node> hostSpeciations = new TreeMap<Double,Node>();
        for (Node node : hostTree.getInternalNodes())
            hostSpeciations.put(node.getHeight(), node);
        
        double[][][] matrices = new double[hostTree.getLeafNodeCount()][][];
        
        for (int i = 0; i < matrices.length; ++i) {
            matrices[i] = constructRateMatrix(i+1);
        }
        
        double density = calculateDensity(originHeight, 1,
                embeddedTree.getRoot(), hostSpeciations, matrices);
        assert(!Double.isNaN(Math.log(density)));
        return Math.log(density);
        
    }
    
    protected double calculateDensity(double originHeight, int startState,
            Node embedded, TreeMap<Double,Node> hostSpeciations,
            double[][][] matrices) {
        
        double L = 1.0;
        
        Tree hostTree = hostTreeInput.get();
        Reconciliation reconciliation = reconciliationInput.get();
        
        BranchRateModel branchRateModel = branchRateModelInput.get();
        double rate = branchRateModel.getRateForBranch(embedded);
        
        double embeddedHeight = embedded.getHeight();
        
        Node host = reconciliation.getHost(embedded);
        if (!Utils.lineageExistedAtHeight(host, embeddedHeight))
            return 0.0;
        
        double startHeight = originHeight;
        
        Map<Node,Integer> hostNodes2Bins = new HashMap<Node,Integer>();
        List<Node> startHostLineages = Utils.getLineagesAtHeight(hostTree,
                startHeight, true);
        int hostCount = startHostLineages.size();
        for (int i = 0; i < hostCount; ++i)
            hostNodes2Bins.put(startHostLineages.get(i), i);
        
        NavigableSet<Double> hostSpeciationSet = hostSpeciations.subMap(
                embeddedHeight, false, startHeight, false).descendingKeySet();
        if (hostSpeciationSet.size() > 0)
            assert(!MachineAccuracy.same(startHeight, hostSpeciationSet.first()));
        double[] startDensity = new double[getStateCount(hostCount)];
        startDensity[startState] = 1.0;
        
        assert(hostNodes2Bins.containsKey(host) || hostSpeciationSet.size() > 0);
        int speciatedBin = -1;
        if (hostSpeciationSet.isEmpty())
            speciatedBin = hostNodes2Bins.get(host);
        double[] endDensity;
        
        for (double speciationHeight : hostSpeciationSet) {
            
            double[][] matrix = matrices[hostCount-1];
            double t = (startHeight - speciationHeight) * rate;
            endDensity = AMH11.expmv(t, matrix, startDensity);
            startHeight = speciationHeight;
            speciatedBin = hostNodes2Bins.get(
                    hostSpeciations.get(speciationHeight));
            hostNodes2Bins.clear();
            startHostLineages = Utils.getLineagesAtHeight(hostTree, startHeight,
                    true);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            ++hostCount;
            assert(startHostLineages.size() == hostCount);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            int[] map = mapNewStatesToOld(hostCount, speciatedBin);
            startDensity = new double[getStateCount(hostCount)];
            for (int i = 0; i < map.length; ++i) {
                double density =
                        map[i] != -1 ? endDensity[i] : 0;
                        startDensity[i] = density;
            }
            
        }

        assert(startHostLineages.size() == hostCount);
        double[][] matrix = matrices[hostCount-1];
        double t = (startHeight - embeddedHeight) * rate;
        endDensity = AMH11.expmv(t, matrix, startDensity);
                
        Node embeddedLeft = embedded.getLeft();
        Node embeddedRight = embedded.getRight();
        int[] state = new int[hostCount];

        if (embedded.isLeaf()) {
            
            int hostBin = hostNodes2Bins.get(host);
//            if (hostSpeciationSet.size() > 0)
//                state = new int[hostCount-1];
            assert(hostBin < state.length);
            state[hostBin] = 1;
//            assert(endDensity.get(compressState(state)) > 0);
            return endDensity[compressState(state)];
            
        } else if (MachineAccuracy.same(embeddedHeight, host.getHeight())) {

            ++hostCount;
            state = new int[hostCount];
            hostNodes2Bins.clear();
            startHostLineages =
                    Utils.getLineagesAtHeight(hostTree, embeddedHeight, true);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            
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
            
            int[] map = mapStatesToOneMore(hostCount-1, speciatedBin);
            int[] map1 = mapNewStatesToOld(hostCount, speciatedBin);
            
            int stateCount = getStateCount(hostCount);
            startDensity =
                    new double[stateCount];
            for (int i = 0; i < stateCount; ++i) {
                double density = map1[i] != -1 && map[map1[i]] != -1 ?
                        endDensity[map[map1[i]]] : 0;
                startDensity[i] = density;
            }
            
            L *= calculateDensity(embeddedHeight, startDensity, hostSpeciations,
                    rate, matrices);
//            assert(L > 0);
            assert(!Double.isNaN(L) && !Double.isInfinite(L));

        } else {
            
            double lambda = duplicationRateParameterInput.get().getValue();
            double tau = hostSwitchRateParameterInput.get().getValue();
            double lambdaPtau = lambda + tau;
            if (lambdaPtau == 0) return 0.0;
            double pDuplication = lambda / lambdaPtau;
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
            
            if (hostCount > 1)
                sum1 += (1 - pDuplication) * sum2 / (hostCount - 1);
            
            state[hostBin] = 1;
            double sum = sum1 * calculateDensity(embeddedHeight,
                    compressState(state), embeddedRight, hostSpeciations,
                    matrices);
            state[hostBin] = 0;

            integrated = embeddedRight;
            sum1 = 0.0;
            for (int i = 0; i < hostCount; ++i) {
                if (i != hostBin) {
                    state[i] = 1;
                    sum1 += calculateDensity(embeddedHeight,
                            compressState(state), integrated, hostSpeciations,
                            matrices);
                    state[i] = 0;
                }
            }

            if (hostCount > 1)
                sum1 *= (1 - pDuplication) * (hostCount - 1);

            state[hostBin] = 1;
            sum += sum1 * calculateDensity(embeddedHeight, compressState(state),
                    embeddedLeft, hostSpeciations, matrices);
            state[hostBin] = 0;
                        
            L *= sum;
            
            int[] map = mapStatesToOneMore(hostCount, hostBin);
            int stateCount = getStateCount(hostCount);
            startDensity = new double[stateCount];
            for (int i = 0; i < stateCount; ++i) {
                double density = map[i] != -1 ? endDensity[map[i]] : 0;
                startDensity[i] = density;
            }
            
            L *= calculateDensity(embeddedHeight, startDensity, hostSpeciations,
                    rate, matrices);
//            assert(L > 0);
            assert(!Double.isNaN(L) && !Double.isInfinite(L));
        }
        
//        assert(L > 0);
        assert(!Double.isNaN(L) && !Double.isInfinite(L));
        return L;
        
    }
    
    protected double calculateDensity(double startHeight,
            double[] startDensity, TreeMap<Double,Node> hostSpeciations,
            double rate, double[][][] matrices) {
        
        Tree hostTree = hostTreeInput.get();
                        
        Map<Node,Integer> hostNodes2Bins = new HashMap<Node,Integer>();
        List<Node> startHostLineages = Utils.getLineagesAtHeight(hostTree,
                startHeight, true);
        int hostCount = startHostLineages.size();
        for (int i = 0; i < hostCount; ++i)
            hostNodes2Bins.put(startHostLineages.get(i), i);
        
        NavigableSet<Double> hostSpeciationSet = hostSpeciations.subMap(0.0,
                false, startHeight, false).descendingKeySet();
                
        int speciatedBin;
        double[] endDensity = new double[hostCount];
        
        for (double speciationHeight : hostSpeciationSet) {
            
            double[][] matrix = matrices[hostCount-1];
            double t = (startHeight - speciationHeight) * rate;
            endDensity = AMH11.expmv(t, matrix, startDensity);
            startHeight = speciationHeight;
            speciatedBin = hostNodes2Bins.get(
                    hostSpeciations.get(speciationHeight));
            hostNodes2Bins.clear();
            startHostLineages = Utils.getLineagesAtHeight(hostTree, startHeight,
                    true);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            ++hostCount;
            assert(startHostLineages.size() == hostCount);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            int[] map = mapNewStatesToOld(hostCount, speciatedBin);
            startDensity = new double[getStateCount(hostCount)];
            for (int i = 0; i < map.length; ++i) {
                double density =
                        map[i] != -1 ? endDensity[map[i]] : 0;
                        startDensity[i] = density;
            }
            
        }

        assert(startHostLineages.size() == hostCount);

        double[][] matrix = matrices[hostCount-1];
        double t = startHeight * rate;
        endDensity = AMH11.expmv(t, matrix, startDensity);
        
//        assert(endDensity.get(0) > 0);
        return endDensity[0];
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
    
    protected double[][] constructRateMatrix(int hostCount) {
        
        double lambda = duplicationRateParameterInput.get().getValue();
        double tau = hostSwitchRateParameterInput.get().getValue();
        double mu = lossRateParameterInput.get().getValue();
        
        int stateCount = getStateCount(hostCount);
        int[][] decompressedStates = new int[stateCount][];
        for (int i = 0; i < stateCount; ++i)
            decompressedStates[i] = decompressState(i, hostCount);
        
        // A sparse matrix would be ideal but dense is substantially faster
        double[][] matrix = new double[stateCount][stateCount];
        
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
                        if (hostCount > 1)
                            rate += (Utils.sum(state) - state[eventLocus]) /
                                    (hostCount - 1) * tau;
                    }
                    
                    eventLocus = locateBirthEvent(decompressedStates[j], state);
                    if (eventLocus != -1)
                        rate += state[eventLocus] * mu;
                    
                }
                
                if (Math.abs(rate) > 0.0) matrix[j][i] = rate;
                
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
        int[] map = new int[newStateCount];
        int[] oldState = new int[hostCount-1];
        for (int state = 0; state < newStateCount; ++state) {
            boolean valid = true;
            decompressState(state, decompressedState);
            for (int i = 0; i < hostCount; ++i) {
                int j = i > speciatedBin ? i-1 : i;
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
        assert(hostCount > 0);
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
 
    protected final int[] mapStatesToOneMore(int hostCount, int speciatedBin) {
        assert(hostCount > 0);
        int newStateCount = getStateCount(hostCount);
        int[] decompressedState = new int[hostCount];
        int[] map = new int[getStateCount(hostCount)];
        for (int state = 0; state < newStateCount; ++state) {
            decompressState(state, decompressedState);
            if (decompressedState[speciatedBin] >= threshold) {
                map[state] = -1;
            } else {
                decompressedState[speciatedBin] += 1;
                map[state] = compressState(decompressedState);
            }
        }
        return map;
    }
    
    protected final int[] mapStatesToNone(int hostCount, int speciatedBin) {
        assert(hostCount > 0);
        int newStateCount = getStateCount(hostCount);
        int[] decompressedState = new int[hostCount];
        int[] map = new int[getStateCount(hostCount)];
        for (int state = 0; state < newStateCount; ++state) {
            decompressState(state, decompressedState);
            decompressedState[speciatedBin] = 0;
            map[state] = compressState(decompressedState);
        }
        return map;       
    }
    
    @Override
    protected double getOverallRate(Node embedded, Node host) {
        
        BranchRateModel branchRateModel = branchRateModelInput.get();
        return branchRateModel.getRateForBranch(embedded);
        
    }
    
}
