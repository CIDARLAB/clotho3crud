/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.JsonToken;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.JsonDeserializer;
import java.io.IOException;
import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.util.ByteSource;
import org.clothocad.core.schema.Constraint;

/**
 *
 * @author spaige
 */
public class SimpleHashDeserializer extends JsonDeserializer<SimpleHash> {

    @Override
    public SimpleHash deserialize(JsonParser jp, DeserializationContext ctxt) throws IOException, JsonProcessingException {
        Integer iterations = null;
        String algorithmName = null;
        byte[] bytes = null;
        ByteSource salt = null;
        JsonToken t = jp.getCurrentToken();
        while (t != JsonToken.END_OBJECT) {
            if (t == JsonToken.FIELD_NAME) {
                String fieldName = jp.getCurrentName();
                t = jp.nextToken();
                switch (fieldName) {
                    case "algorithmName":
                        if (t != JsonToken.VALUE_STRING) {
                            throw ctxt.wrongTokenException(jp, JsonToken.VALUE_STRING, "Expected string specifying algorithm name");
                        }

                        algorithmName = jp.getValueAsString();
                        break;
                    case "bytes":
                        if (t != JsonToken.VALUE_EMBEDDED_OBJECT) {
                            throw ctxt.wrongTokenException(jp, JsonToken.START_OBJECT, "Expected byte source for hashed value");
                        }

                        bytes = (byte[]) jp.getEmbeddedObject();
                        break;
                    case "iterations":
                        if (t != JsonToken.VALUE_NUMBER_INT) {
                            throw ctxt.wrongTokenException(jp, JsonToken.VALUE_NUMBER_INT, "Expected integer number of iterations");
                        }
                        iterations = jp.getIntValue();
                        break;
                    case "salt":
                        if (t != JsonToken.START_OBJECT) {
                            throw ctxt.wrongTokenException(jp, JsonToken.VALUE_EMBEDDED_OBJECT, "Expected byte source for salt");
                        }
                        t = jp.nextToken();
                        if (t != JsonToken.FIELD_NAME && jp.getCurrentName() == "bytes") {
                            throw ctxt.wrongTokenException(jp, JsonToken.FIELD_NAME, "Expected bytes field");
                        }
                        t = jp.nextToken();
                        if (t != JsonToken.VALUE_EMBEDDED_OBJECT) {
                            throw ctxt.wrongTokenException(jp, JsonToken.VALUE_EMBEDDED_OBJECT, "Expected byte source value for salt");
                        }

                        salt = ByteSource.Util.bytes((byte[]) jp.getEmbeddedObject());
                        while (jp.nextToken() != JsonToken.END_OBJECT){
                            
                        }
                        break;
                    default:
                        throw ctxt.weirdStringException(jp.getCurrentName(), Constraint.class, "Must be \"algorithmName\", \"bytes\", \"iterations\", or \"salt\"");

                }
            } else {
                throw ctxt.wrongTokenException(jp, JsonToken.FIELD_NAME, "Expected field name named \"algorithmName\", \"bytes\", \"iterations\", or \"salt\"");
            }
            t = jp.nextToken();
        }
        if (algorithmName == null || iterations == null || bytes == null || salt == null) {
            throw ctxt.instantiationException(SimpleHash.class, "Missing field data");
        }


        SimpleHash hash = new SimpleHash(algorithmName);
        hash.setBytes(bytes);
        hash.setSalt(salt);
        hash.setIterations(iterations);
        return hash;

    }
}
