/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonTypeInfo;
import com.fasterxml.jackson.annotation.JsonValue;
import lombok.EqualsAndHashCode;

/**
 *
 * @author spaige
 */
@EqualsAndHashCode
@JsonTypeInfo(use = JsonTypeInfo.Id.NONE)
public class ObjectId {
    
    @JsonCreator
    public ObjectId(String value){
        this.value = value;
    }
    
    public ObjectId(Object value){
        if (value == null) throw new IllegalArgumentException();
        this.value = value.toString();
    }
    
    public ObjectId(){
        this.value = new org.bson.types.ObjectId().toString();
    }
    
    @JsonValue
    public String getValue(){
        return value;
    }

    @Override
    public String toString() {
        return value;
    }
    
    private final String value;
}
