/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.JsonToken;
import com.fasterxml.jackson.databind.BeanProperty;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.JsonDeserializer;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.deser.ContextualDeserializer;
import com.fasterxml.jackson.databind.jsontype.TypeDeserializer;
import java.io.IOException;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
//TODO: figure out what to do when a reference can't be resolved
public class ReferenceDeserializer extends JsonDeserializer<ObjBase> implements ContextualDeserializer {

    BeanProperty property;

    @Override
    public ObjBase deserialize(JsonParser jp, DeserializationContext ctxt) throws IOException, JsonProcessingException {
        //get id token
        JsonToken t = jp.getCurrentToken();
        if (t.isScalarValue()){
        //get injected value from context
            return (ObjBase) ctxt.findInjectableValue(jp.getValueAsString(), property, null);
        }
        else throw new IllegalArgumentException("Cannot read objectid from "+t);
    }

    @Override
    public Object deserializeWithType(JsonParser jp, DeserializationContext ctxt, TypeDeserializer typeDeserializer) throws IOException, JsonProcessingException {
        return deserialize(jp, ctxt);
    }

    @Override
    public JsonDeserializer<?> createContextual(DeserializationContext ctxt, BeanProperty property) throws JsonMappingException {
        this.property = property;
        return this;
    }
    
}
