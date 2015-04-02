/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.persistence.Persistor;

/**
 *
 * @author spaige
 */
public class MessageOptions {
    public MessageOptions(Map<MessageOption, Object> options){
        if (options != null) this.options = options;
        else this.options = new HashMap<>();
    }
    
    public MessageOptions(){
        this.options = new HashMap<>();
    }
    
    private final Map<MessageOption, Object> options;
    
    public boolean isMute(){
        Boolean value = convert(MessageOption.mute, Boolean.class);
        if (value != null){
            return value;
        }
        return false;
    };
    
    public Set<String> getPropertiesFilter(){
        List<String> value = convert(MessageOption.filter, List.class);
        if (value != null){
            
            Set<String> filter = new HashSet<>();
            filter.addAll(value);
            return filter;
        }
        return new HashSet<>();
    }
    
    public Detail getDetail(){
        String value = convert(MessageOption.detail, String.class);
        if (value != null){
            try {
                return Detail.valueOf(value);
            } catch (IllegalArgumentException e) {
                //pass through and return the default
            }
        }
        return Detail.NORMAL;
    }
   
    public int getMaxResults(){
        Integer value = convert(MessageOption.maxResults, Integer.class);
        if (value != null){
            return value;
        }
        return Persistor.SEARCH_MAX;
    }
    
    
    private <T> T convert(MessageOption key, Class<T> c){
        if (options.containsKey(key)){
            Object value = options.get(key);
            if (c.isInstance(value)) return (T) value;
        }
        return null;
    }
    
    public static enum Detail {
        NORMAL,
        ID_ONLY
    }
}
