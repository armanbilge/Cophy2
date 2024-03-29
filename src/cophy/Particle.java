/**
 * Particle.java
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
import java.util.Map;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class Particle<T> implements Cloneable {

    protected T value;
    protected double weight;
    
    public Particle() {}
    
    public Particle(T value) {
        this(value, 1.0);
    }
    
    public Particle(T value, double weight) {
        this.value = value;
        this.weight = weight;
    }
    
    public T getValue() {
        return value;
    }
    
    public void setValue(T value) {
        this.value = value;
    }
    
    public double getWeight() {
        return weight;
    }
    
    public void resetWeight() {
        weight = 1.0;
    }
    
    public void muliWeight(double value) {
        weight *= value;
    }
    
    @Override
    public Particle<T> clone() {
        return new Particle<T>(value, weight);
    }
    
    public static class TreeParticle extends Particle<Tree> {
        
        protected Map<Node,Node> nodeMap;
        
        {
            nodeMap = new HashMap<Node,Node>();
        }
        
        public TreeParticle() {
            super();
        }
        
        public TreeParticle(Tree value) {
            super(value);
        }
        
        public TreeParticle(Tree value, double weight) {
            super(value, weight);
        }
        
        public Map<Node,Node> getNodeMap() {
            return nodeMap;
        }
        
        @Override
        public TreeParticle clone() {
            TreeParticle clone = (TreeParticle) super.clone();
            clone.nodeMap = new HashMap<Node,Node>(nodeMap);
            return clone;
        }
        
    }
    
}
