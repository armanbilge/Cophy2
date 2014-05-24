/**
 * CarryingCapacityCophylogenyModel.java
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

import org.apache.commons.math.util.MathUtils;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class CarryingCapacityCophylogenyModel extends EmbeddedTreeDistribution {
    
    protected int carryingCapacity;
    
    /**
     * 
     */
    public CarryingCapacityCophylogenyModel() {
        // TODO Auto-generated constructor stub
    }

    protected DoubleMatrix2D constructRateMatrix(int hostCount, double lambda,
            double tau, double mu) {
        
        int stateCount = MathUtils.pow(carryingCapacity + 1, hostCount);
        
        DoubleFactory2D.sparse.make(stateCount, stateCount);
        
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                
                
                
            }
        }
        
        return null;
    }
    
    protected final int[] decomposeState(int state, int hostCount) {
        
        int base = carryingCapacity + 1;
        int[] decomposedState = new int[hostCount];
        for (int i = 0; i < hostCount; ++i) {
            decomposedState[i] = state % base;
            state /= base;
        }
        return decomposedState;
        
    }
    
    protected final boolean createdBy1BirthEvent(int[] s0, int[] s1) {
        
        boolean eventCounted = false;
        
        for (int i = 0; i < s0.length; ++i) {
            int change = s1[i] - s0[i];
            if (change > 1 || change < 0) return false;
            if (change == 1 && eventCounted) return false;
            else eventCounted = true;
        }
        
        return true;
        
    }
    
}
