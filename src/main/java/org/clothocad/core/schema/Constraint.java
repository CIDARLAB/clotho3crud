/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonTypeInfo;
import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.BeanProperty;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.JsonDeserializer;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.deser.ContextualDeserializer;
import com.fasterxml.jackson.databind.jsontype.TypeDeserializer;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author spaige
 */
public class Constraint {
    public final Class constraintType;
    @JsonTypeInfo(use = JsonTypeInfo.Id.NONE)
    public final Map<String, Object> values;
    
    public Constraint(@JsonProperty("constraintType") Class constraintType, @JsonProperty("values") Map<String,Object> values){
        this.values = values;
        this.constraintType = constraintType;
    }
    
    public Constraint(Class annotation, Object... args){
        this.constraintType = annotation;
        values = new HashMap<>();
        if (args.length % 2 != 0){
            throw new IllegalArgumentException("Needs value for each value name");
        }
        for (int i=0; i< args.length; i=i+2){
            if (!(args[i] instanceof String)){
                throw new IllegalArgumentException("Expected string value name");
            }
            values.put((String) args[i], args[i+1]);
        }
    }

}
