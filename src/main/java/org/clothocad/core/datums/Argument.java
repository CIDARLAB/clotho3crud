package org.clothocad.core.datums;

import org.clothocad.core.persistence.IdUtils;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.JsonToken;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.JsonDeserializer;
import com.fasterxml.jackson.databind.JsonSerializer;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.annotation.JsonSerialize;
import com.fasterxml.jackson.databind.jsontype.TypeSerializer;
import com.fasterxml.jackson.databind.util.TokenBuffer;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

@NoArgsConstructor
@Slf4j
public class Argument {

    @Getter
    private String name;
    @Getter
    //XXX: Should probably be Type instead
    @JsonSerialize(using = ArgTypeSerializer.class)
    @JsonDeserialize(using = ArgTypeDeserializer.class)
    private Class type;

    public Argument(String name, Class type) {
        this.name = name;
        this.type = type;
    }

    public static String jsonifyFieldType(Class c) {

        if (ObjBase.class.isAssignableFrom(c)) {
            //we should write the schema name
            return c.getName();
        }
        if (ObjectId.class.isAssignableFrom(c)) {
            return "id";
        }
        if (Date.class.isAssignableFrom(c)) {
            return "date";
        }
        if (String.class.isAssignableFrom(c) || c.equals(char.class)) {
            return "string";
        }
        if (Boolean.class.isAssignableFrom(c) || c.equals(boolean.class)) {
            return "boolean";
        }
        if (Number.class.isAssignableFrom(c)
                || c.equals(byte.class)
                || c.equals(short.class)
                || c.equals(int.class)
                || c.equals(long.class)
                || c.equals(float.class)
                || c.equals(double.class)) {
            return "number"; // String.format("number(%s)", c.getSimpleName());
        }
        if (c.isArray() || Collection.class.isAssignableFrom(c)) {
            //todo: parameterize array types;
            return "array";
        }
        if (Map.class.isAssignableFrom(c)) {
            return "object";
        }
        log.warn("Unable to jsonify field type {}", c.getName());
        return "object";
    }

    protected final static Map<String, Class> classMap;

    static {
        //todo: make immutable
        classMap = new HashMap<>();
        classMap.put("string", String.class);
        classMap.put("boolean", Boolean.class);
        classMap.put("number", Number.class);
        classMap.put("array", List.class);
        classMap.put("id", ObjectId.class);
        classMap.put("date", Date.class);
        classMap.put("object", Map.class);
    }

    public static Class decodeFieldType(String s) {
        Class c = classMap.get(s.toLowerCase());
        if (c != null) {
            return c;
        }
        ObjectId id = new ObjectId(s);
        try {
            return IdUtils.getClass(id);
        } catch (ClassNotFoundException ex) {
            throw new RuntimeException("Could not find schema: " + s, ex);
        }
    }

    public static class ArgTypeSerializer extends JsonSerializer<Class> {

        @Override
        public void serialize(Class value, JsonGenerator jgen, SerializerProvider provider) throws IOException, JsonProcessingException {
            jgen.writeString(jsonifyFieldType(value));
        }

        @Override
        public void serializeWithType(Class value, JsonGenerator jgen, SerializerProvider provider, TypeSerializer typeSer) throws IOException, JsonProcessingException {
            serialize(value, jgen, provider);
        }
    }

    public static class ArgTypeDeserializer extends JsonDeserializer<Class> {

        @Override
        public Class deserialize(JsonParser jp, DeserializationContext ctxt) throws IOException, JsonProcessingException {
            JsonToken t = jp.getCurrentToken();
            if (t.isScalarValue()) {
                return decodeFieldType(jp.getValueAsString());
            } else {
                throw new IllegalArgumentException("Cannot read argument type from " + t);
            }
        }
    }
}
