/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JsonSerializer;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.jsontype.TypeSerializer;
import java.io.IOException;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public class ReferenceSerializer extends JsonSerializer<ObjBase>{

    @Override
    public void serialize(ObjBase value, JsonGenerator jgen, SerializerProvider provider) throws IOException, JsonProcessingException {
        if (value.getId() == null){
            value.setId(new ObjectId());
        }
        ObjectId id = value.getId();
        jgen.writeString(id.toString());
    }   

    @Override
    public void serializeWithType(ObjBase value, JsonGenerator jgen, SerializerProvider provider, TypeSerializer typeSer) throws IOException, JsonProcessingException {
        //serialization w/ type is irrelevant - type is encoded in referenced document
        serialize(value, jgen, provider);
    }
    
}