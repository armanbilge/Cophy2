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

import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.State;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public abstract class EmbeddedTreeDistribution extends Distribution {

    /**
     * 
     */
    public EmbeddedTreeDistribution() {
        // TODO Auto-generated constructor stub
    }

    /* (non-Javadoc)
     * @see beast.core.Distribution#getArguments()
     */
    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    /* (non-Javadoc)
     * @see beast.core.Distribution#getConditions()
     */
    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    /* (non-Javadoc)
     * @see beast.core.Distribution#sample(beast.core.State, java.util.Random)
     */
    @Override
    public void sample(State state, Random random) {
        // TODO Auto-generated method stub

    }

}
