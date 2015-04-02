/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.annotation.JsonAnySetter;

import java.util.HashMap;
import java.util.Map;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public class JSONFilter extends ObjBase{
    
    //@JsonSerialize(using=FilterSerializer.class)
    private Map<String,Object> _fields = new HashMap<>();
    
    @JsonAnySetter
    public void stashField(String name, Object o){
        if (name.equals("_id")) name = ID;
        
        
        _fields.put(name, o);
    }
    
    @Override
    public ObjectId getId(){
        if (super.getId() != null) return super.getId();
        Object id = _fields.get(ID);
        if (id == null) return null;
        return new ObjectId(id);
    }
    
    
    /*
    public static class FilterSerializer extends JsonSerializer<Object> {

        @Override
        public void serialize(Object value, JsonGenerator jgen, SerializerProvider provider) throws IOException, JsonProcessingException {
            provider.g
        }

        @Override
        public void serializeWithType(Object value, JsonGenerator jgen, SerializerProvider provider, TypeSerializer typeSer) throws IOException, JsonProcessingException {
            serialize(value, jgen, provider);
        }
        
        
    }*/
}
