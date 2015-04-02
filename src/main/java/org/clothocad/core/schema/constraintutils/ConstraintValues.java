/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema.constraintutils;

import java.util.Set;

/**
 *
 * @author spaige
 */
public interface ConstraintValues {
    public Object get(String key);
    
    public void put(String key, Object value);
    
    public Set<String> keySet();   
}
