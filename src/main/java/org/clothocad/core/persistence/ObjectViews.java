/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import java.util.HashMap;
import java.util.Map;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
public class ObjectViews {
    public <T extends ObjBase> T get(Class<T> type){
        return (T) map.get(type);
    }
    
    public <T extends ObjBase> void put(Class<T> type, T instance){
        map.put(type, instance);
    }
    
    public boolean has(Class<? extends ObjBase> type){
        return map.containsKey(type);
    }
    
    Map<Class<? extends ObjBase>, ObjBase> map = new HashMap();
}
