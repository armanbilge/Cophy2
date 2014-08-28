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

import static org.apache.commons.math.special.Gamma.logGamma;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeMap;

import org.apache.commons.math.util.MathUtils;

import amh11.AMH11;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cophy.Reconciliation;
import cophy.Util;

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
        
        DoubleMatrix2D[] matrices =
                new DoubleMatrix2D[hostTree.getLeafNodeCount()];
        
        for (int i = 0; i < matrices.length; ++i) {
            matrices[i] = constructRateMatrix(i+1);
        }
        
        double density = calculateDensity(originHeight, 1,
                embeddedTree.getRoot(), hostSpeciations, matrices);
        assert(!Double.isNaN(Math.log(density)));
        int taxonCount = embeddedTree.getLeafNodeCount();
        double logL = Math.log(density)
                + (taxonCount - 1) * Math.log(2.0)
                - logGamma(taxonCount + 1) - logGamma(taxonCount);
        return logL;
    }
    
    protected double calculateDensity(double originHeight, int startState,
            Node embedded, TreeMap<Double,Node> hostSpeciations,
            DoubleMatrix2D[] matrices) {
        
        double L = 1.0;
        
        Tree hostTree = hostTreeInput.get();
        Reconciliation reconciliation = reconciliationInput.get();
        
        double rate = rateParameterInput.get().getValue();
        
        double embeddedHeight = embedded.getHeight();
        
        Node host = reconciliation.getHost(embedded);
        if (embeddedHeight == originHeight
                || !Util.lineageExistedAtHeight(host, embeddedHeight))
            return 0.0;
                
        double startHeight = originHeight;
        
        Map<Node,Integer> hostNodes2Bins = new HashMap<Node,Integer>();
        List<Node> startHostLineages = Util.getLineagesAtHeight(hostTree,
                startHeight, true);
        int hostCount = startHostLineages.size();
        for (int i = 0; i < hostCount; ++i)
            hostNodes2Bins.put(startHostLineages.get(i), i);
        
        NavigableSet<Double> hostSpeciationSet = hostSpeciations.subMap(
                embeddedHeight, false, startHeight, false).descendingKeySet();
        if (hostSpeciationSet.size() > 0)
            assert(startHeight != hostSpeciationSet.first());
        DoubleMatrix1D startDensity =
                DoubleFactory1D.dense.make(getStateCount(hostCount));
        startDensity.setQuick(startState, 1.0);
        
        assert(hostNodes2Bins.containsKey(host) || hostSpeciationSet.size() > 0);
        int speciatedBin = -1;
        if (hostSpeciationSet.isEmpty())
            speciatedBin = hostNodes2Bins.get(host);
        DoubleMatrix1D endDensity;
        
        for (double speciationHeight : hostSpeciationSet) {
            
            DoubleMatrix2D matrix = matrices[hostCount-1];
            double t = (startHeight - speciationHeight) * rate;
            assert(matrix.columns() == startDensity.size());
            endDensity = AMH11.expmv(t, matrix, startDensity);
            startHeight = speciationHeight;
            speciatedBin = hostNodes2Bins.get(
                    hostSpeciations.get(speciationHeight));
            hostNodes2Bins.clear();
            startHostLineages = Util.getLineagesAtHeight(hostTree, startHeight,
                    true);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            ++hostCount;
            assert(startHostLineages.size() == hostCount);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            int[] map = mapNewStatesToOld(hostCount, speciatedBin);
            startDensity = DoubleFactory1D.dense.make(getStateCount(hostCount));
            for (int i = 0; i < map.length; ++i) {
                double density =
                        map[i] != -1 ? endDensity.getQuick(map[i]) : 0;
                        startDensity.setQuick(i, density);
            }
            
        }

        assert(startHostLineages.size() == hostCount);
        DoubleMatrix2D matrix = matrices[hostCount-1];
        double t = (startHeight - embeddedHeight) * rate;
        assert(matrix.columns() == startDensity.size());
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
            assert(compressState(state) < endDensity.size());
//            assert(endDensity.getQuick(compressState(state)) > 0);
            assert(!Double.isNaN(endDensity.getQuick(compressState(state)))
                    && !Double.isInfinite(endDensity.getQuick(compressState(state))));
            return endDensity.getQuick(compressState(state));
            
        } else if (embeddedHeight == host.getHeight()) {

            ++hostCount;
            state = new int[hostCount];
            hostNodes2Bins.clear();
            startHostLineages =
                    Util.getLineagesAtHeight(hostTree, embeddedHeight, true);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            
            int hostBin = hostNodes2Bins.get(host.getLeft());
            state[hostBin] = 1;
            int compressedState = compressState(state);
            double p1 = calculateDensity(embeddedHeight, compressedState,
                    embeddedLeft, hostSpeciations, matrices);
            double p2 = calculateDensity(embeddedHeight, compressedState,
                    embeddedRight, hostSpeciations, matrices);
            state[hostBin] = 0;

            hostBin = hostNodes2Bins.get(host.getRight());
            state[hostBin] = 1;
            compressedState = compressState(state);
            p1 *= calculateDensity(embeddedHeight, compressedState,
                    embeddedRight, hostSpeciations, matrices);
            p2 *= calculateDensity(embeddedHeight, compressedState,
                    embeddedLeft, hostSpeciations, matrices); 
            state[hostBin] = 0;
            
            L *= p1 + p2;
            
            int[] map = mapStatesToOneMore(hostCount-1, speciatedBin);
            int[] map1 = mapNewStatesToOld(hostCount, speciatedBin);
            
            int stateCount = getStateCount(hostCount);
            startDensity = DoubleFactory1D.dense.make(stateCount);
            for (int i = 0; i < stateCount; ++i) {
                double density = map1[i] != -1 && map[map1[i]] != -1 ?
                        endDensity.getQuick(map[map1[i]]) : 0;
                startDensity.setQuick(i, density);
            }
            
            L *= calculateDensity(embeddedHeight, startDensity, hostSpeciations,
                    rate, matrices);
//            assert(L > 0);
            assert(!Double.isNaN(L) && !Double.isInfinite(L));

        } else {
            
            double lambda = duplicationRateParameterInput.get().getValue();
            double tau = hostSwitchRateParameterInput.get().getValue();
            int hostBin = hostNodes2Bins.get(host);
            
            for (int i = 0; i < getStateCount(hostCount); ++i)
                endDensity.setQuick(i, endDensity.getQuick(i)
                        * decompressState(i, hostCount)[hostBin]);

            state[hostBin] = 1;
            int compressedState = compressState(state);
            double pSameHostLeft = calculateDensity(embeddedHeight,
                    compressedState, embeddedLeft, hostSpeciations, matrices);
            double pSameHostRight = calculateDensity(embeddedHeight,
                    compressedState, embeddedRight, hostSpeciations, matrices);
            state[hostBin] = 0;
            
            double pIntegratedLeft = 0.0;
            double pIntegratedRight = 0.0;
            if (hostCount > 1) {
                for (int i = 0; i < hostCount; ++i) {
                    if (i != hostBin) {
                        state[i] = 1;
                        compressedState = compressState(state);
                        pIntegratedLeft += calculateDensity(embeddedHeight,
                                compressedState, embeddedLeft,
                                hostSpeciations, matrices);
                        pIntegratedRight += calculateDensity(embeddedHeight,
                                compressedState, embeddedRight,
                                hostSpeciations, matrices);
                        state[i] = 0;
                    }
                }
                pIntegratedLeft /= hostCount - 1;
                pIntegratedRight /= hostCount - 1;
            }
            
            L *= 2 * lambda * pSameHostLeft * pSameHostRight +
                    tau * (pSameHostLeft * pIntegratedRight
                            + pSameHostRight * pIntegratedLeft);

            int[] map = mapStatesToOneMore(hostCount, hostBin);
            int stateCount = getStateCount(hostCount);
            startDensity = DoubleFactory1D.dense.make(stateCount);
            for (int i = 0; i < stateCount; ++i) {
                double density = map[i] != -1 ? endDensity.getQuick(map[i]) : 0;
                startDensity.setQuick(i, density);
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
            DoubleMatrix1D startDensity, TreeMap<Double,Node> hostSpeciations,
            double rate, DoubleMatrix2D[] matrices) {
        
        Tree hostTree = hostTreeInput.get();
                        
        Map<Node,Integer> hostNodes2Bins = new HashMap<Node,Integer>();
        List<Node> startHostLineages = Util.getLineagesAtHeight(hostTree,
                startHeight, true);
        int hostCount = startHostLineages.size();
        for (int i = 0; i < hostCount; ++i)
            hostNodes2Bins.put(startHostLineages.get(i), i);
        
        NavigableSet<Double> hostSpeciationSet = hostSpeciations.subMap(0.0,
                false, startHeight, false).descendingKeySet();
                
        int speciatedBin;
        DoubleMatrix1D endDensity;
        
        for (double speciationHeight : hostSpeciationSet) {
            
            DoubleMatrix2D matrix = matrices[hostCount-1];
            double t = (startHeight - speciationHeight) * rate;
            assert(matrix.columns() == startDensity.size());
            endDensity = AMH11.expmv(t, matrix, startDensity);
            startHeight = speciationHeight;
            speciatedBin = hostNodes2Bins.get(
                    hostSpeciations.get(speciationHeight));
            hostNodes2Bins.clear();
            startHostLineages = Util.getLineagesAtHeight(hostTree, startHeight,
                    true);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            ++hostCount;
            assert(startHostLineages.size() == hostCount);
            for (int i = 0; i < hostCount; ++i)
                hostNodes2Bins.put(startHostLineages.get(i), i);
            int[] map = mapNewStatesToOld(hostCount, speciatedBin);
            startDensity = DoubleFactory1D.dense.make(getStateCount(hostCount));
            for (int i = 0; i < map.length; ++i) {
                double density =
                        map[i] != -1 ? endDensity.getQuick(map[i]) : 0;
                        startDensity.setQuick(i, density);
            }
            
        }

        assert(startHostLineages.size() == hostCount);

        DoubleMatrix2D matrix = matrices[hostCount-1];
        double t = startHeight * rate;
        assert(matrix.columns() == startDensity.size());
        endDensity = AMH11.expmv(t, matrix, startDensity);
        
//        assert(endDensity.getQuick(0) > 0);
        assert(!Double.isNaN(endDensity.getQuick(0)) &&
                !Double.isInfinite(endDensity.getQuick(0)));
        return endDensity.getQuick(0);
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
        
        // A sparse matrix would be ideal but dense is substantially faster
        DoubleMatrix2D matrix =
                DoubleFactory2D.dense.make(stateCount, stateCount);
        
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                
                double rate = 0.0;
                int[] state = decompressedStates[i];
                
                if (i != j) {
                    
                    int eventLocus = locateBirthEvent(state,
                            decompressedStates[j]);
                    if (eventLocus != -1) {
                        rate += state[eventLocus] * lambda;
                        if (hostCount > 1)
                            rate += (Util.sum(state) - state[eventLocus]) /
                                    (hostCount - 1) * tau;
                    }
                    
                    eventLocus = locateBirthEvent(decompressedStates[j], state);
                    if (eventLocus != -1)
                        rate += state[eventLocus] * mu;
                    
                }
                
                if (Math.abs(rate) > 0.0) matrix.setQuick(j, i, rate);
                
            }
            
        }
        
        for (int i = 0; i < stateCount; ++i) {
            matrix.set(i, i, -matrix.viewColumn(i).zSum());
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
        
}
