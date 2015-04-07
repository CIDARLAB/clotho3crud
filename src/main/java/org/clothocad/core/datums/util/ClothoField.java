package org.clothocad.core.datums.util;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;
import org.clothocad.core.schema.Access;
import org.clothocad.core.schema.Constraint;
import org.clothocad.core.schema.Schema;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.JsonToken;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.JavaType;
import com.fasterxml.jackson.databind.JsonDeserializer;
import com.fasterxml.jackson.databind.JsonSerializer;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.annotation.JsonSerialize;
import com.fasterxml.jackson.databind.jsontype.TypeDeserializer;
import com.fasterxml.jackson.databind.util.TokenBuffer;

import lombok.Getter;
import lombok.Setter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.GenericArrayType;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.lang.reflect.TypeVariable;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

@Getter
@Setter
public class ClothoField {

    private ClothoField() {}

    public ClothoField(Field field) {

        if (field.getAnnotation(JsonProperty.class) != null) {
            name = field.getAnnotation(JsonProperty.class).value();
        } else {
            name = field.getName();
        }

        type = field.getType();
        reference = field.getAnnotation(Reference.class) != null;
        referenceCollection = field.getAnnotation(ReferenceCollection.class) != null;

        //TODO: access, validate

        //TODO: metadata
        //example
        //description
    }

    private static final Logger logger = LoggerFactory.getLogger(ClothoField.class);

    public ClothoField(String name, Class type, String example, String description, boolean reference, Access access) {
        this.name = name;
        this.type = type;
        this.example = example;
        this.access = access; 
        this.reference = reference;
        this.description = description;
    }

    private String name;

    @Getter
    @JsonSerialize(using = ClothoFieldTypeSerializer.class)
    @JsonDeserialize(using = ClothoFieldTypeDeserializer.class)
    private Class<?> type;
    private Type subtype;
    private String example;   //A string representation/explanation of an expected value
    private Access access;
    private boolean reference;
    private boolean referenceCollection;

    //@JsonSerialize(using=ConstraintsSerializer.class)
    @JsonDeserialize(contentUsing = ConstraintDeserializer.class)
    //@JsonTypeInfo(use = JsonTypeInfo.Id.CLASS, property = "constraintType")
    private Set<Constraint> constraints;

    //metadata
    private String description;

    public String getSetterName() {
        return "set" + capitalize(name);
    }

    public String getGetterName() {
        if (this.type.equals(Boolean.class)) {
            return "is" + capitalize(name);
        }
        return "get" + capitalize(name);
    }

    private static String capitalize(String s) {
        if (s.length() == 0) return s;
        return s.substring(0,1).toUpperCase() + s.substring(1);
    }

    //XXX: needs more specific type info (parameterization) in some cases
    public static class ClothoFieldTypeDeserializer extends JsonDeserializer<Class<?>> {

        @Override
        public Class<?> deserialize(JsonParser jp, DeserializationContext ctxt) throws IOException, JsonProcessingException {
            if (jp.getCurrentToken() != JsonToken.VALUE_STRING) {
                throw ctxt.wrongTokenException(jp, JsonToken.VALUE_STRING, "Expected type name as string");
            }
            String typeName = jp.getValueAsString();
            switch (typeName) {
                case "id":
                    return ObjectId.class;
                case "date":
                    return Date.class;
                case "string":
                    return String.class;
                case "boolean":
                    return Boolean.class;
                case "number":
                    return Number.class;
                case "object":
                    return Map.class;
                case "array":
                    return Collection.class;
                default:
                    try {
                        return Class.forName(typeName, true, Schema.cl);
                    } catch (ClassNotFoundException ex) {
                        throw ctxt.weirdStringException(Class.class, typeName + " is not a schema id or field type");
                    }
            }
        }

        @Override
        public Object deserializeWithType(JsonParser jp, DeserializationContext ctxt, TypeDeserializer typeDeserializer) throws IOException, JsonProcessingException {
            return deserialize(jp,ctxt);
        }

    }

    public static class ClothoFieldTypeSerializer extends JsonSerializer<Class<?>> {

        @Override
        public void serialize(Class<?> value, JsonGenerator jgen, SerializerProvider provider) throws IOException, JsonProcessingException {
            jgen.writeString(jsonifyFieldType(value));
        }

    }

    public static String jsonifyFieldType(Class c) {

        if (ObjBase.class.isAssignableFrom(c)) return c.getCanonicalName();
        if (ObjectId.class.isAssignableFrom(c)) return "id";
        if (Date.class.isAssignableFrom(c)) return "date";
        if (String.class.isAssignableFrom(c) || c.equals(char.class)) return "string";
        if (Boolean.class.isAssignableFrom(c) || c.equals(boolean.class)) return "boolean";
        if (Number.class.isAssignableFrom(c) ||
                c.equals(byte.class) ||
                c.equals(short.class) ||
                c.equals(int.class) ||
                c.equals(long.class) ||
                c.equals(float.class) ||
                c.equals(double.class)) return "number"; //String.format("number(%s)", c.getSimpleName());
        if (c.isArray() || Collection.class.isAssignableFrom(c)) {
           //todo: parameterize array types;
            return "array";
        }
        if (Map.class.isAssignableFrom(c)) return "object";
        logger.warn("Unable to jsonify field type {}", c.getName());
        return "object";
    }
    //Constraints

    //#
    //multipleof
    //maximum
    //exclusivemaximum
    //minimum
    //exclusiveminimum

    //size
    //pattern (regex match)


    //notnull

    public static Type getParameterizedType(Type type) {
        int index = 0;
        if (type instanceof ParameterizedType) {
            ParameterizedType ptype = (ParameterizedType) type;
            if ((ptype.getActualTypeArguments() != null) && (ptype.getActualTypeArguments().length <= index)) {
                return null;
            }
            Type paramType = ptype.getActualTypeArguments()[index];
            if (paramType instanceof GenericArrayType) {
                return ((GenericArrayType) paramType).getGenericComponentType();
            } else {
                if (paramType instanceof ParameterizedType) {
                    return paramType;
                } else {
                    if (paramType instanceof TypeVariable) {
                        // TODO: Figure out what to do... Walk back up the to
                        // the parent class and try to get the variable type
                        // from the T/V/X
//                      throw new MappingException("Generic Typed Class not supported:  <" + ((TypeVariable) paramType).getName() + "> = " + ((TypeVariable) paramType).getBounds()[0]);
                        return paramType;
                    } else if (paramType instanceof Class) {
                        return (Class) paramType;
                    } else {
                        throw new RuntimeException("Unknown type... pretty bad... call for help, wave your hands... yeah!");
                    }
                }
            }
        }
        return null;
    }


    public static class ConstraintDeserializer extends JsonDeserializer<Constraint> {

        @Override
        public Constraint deserialize(JsonParser jp, DeserializationContext ctxt) throws IOException, JsonProcessingException {
            JsonToken t = jp.nextToken();
            Class constraintType = null;
            TokenBuffer tokens = new TokenBuffer(jp.getCodec());
            while (t != JsonToken.END_OBJECT) {
                if (t == JsonToken.FIELD_NAME) {
                    switch (jp.getCurrentName()) {
                        case "constraintType":
                            t = jp.nextToken();
                            if (t != JsonToken.VALUE_STRING) {
                                throw ctxt.wrongTokenException(jp, JsonToken.VALUE_STRING, "Expected string specifying constraint type");
                            }
                            try {
                                //any magic names we want should be handled here
                                constraintType = Schema.cl.loadClass(jp.getValueAsString());
                            } catch (ClassNotFoundException ex) {
                                throw ctxt.weirdStringException(jp.getValueAsString(), Class.class, "Not a valid constraint id.");
                            }
                            break;
                        case "values":
                            t = jp.nextToken();
                            if (t != JsonToken.START_OBJECT) {
                                throw ctxt.wrongTokenException(jp, JsonToken.START_OBJECT, "Expected object specifying constraint parameter values");
                            }
                            tokens.copyCurrentStructure(jp);
                            break;
                        default:
                            throw ctxt.weirdStringException(jp.getCurrentName(), Constraint.class, "Must be \"constraintType\" or \"values\"");
                    }
                } else {
                    throw ctxt.wrongTokenException(jp, JsonToken.FIELD_NAME, "Expected field name named \"constraintType\" or \"values\"");
                }
                t = jp.nextToken();
            }
            //cannot create constraint without a declared type
            if (constraintType == null)  {}
            //use annotation class to find out property types
            //ignore properties not in annotation class
            Map<String,Object> values = new HashMap<>();
            jp = tokens.asParser();
            t = jp.nextToken();

            if (jp.hasCurrentToken()) {
                t = jp.nextToken();
                while (t != JsonToken.END_OBJECT) {
                    if (t != JsonToken.FIELD_NAME) throw ctxt.wrongTokenException(jp, JsonToken.FIELD_NAME, "Constraint parameter values malformed");
                    String name = jp.getCurrentName();
                    jp.nextToken();
                    try {
                        JavaType fieldType = ctxt.getTypeFactory().constructType(constraintType.getMethod(name).getReturnType());
                        Object fieldValue = ctxt.findContextualValueDeserializer(fieldType, null).deserialize(jp, ctxt);
                        values.put(name, fieldValue);
                    } catch (NoSuchMethodException ex) {
                        jp.skipChildren();
                    } catch (SecurityException ex) {
                        throw new RuntimeException("Couldn't access field " + jp.getCurrentName() + " in " + constraintType.getCanonicalName());
                    }
                    t = jp.nextToken();
                }
            }

            return new Constraint(constraintType, values);
        }

        @Override
        public Object deserializeWithType(JsonParser jp, DeserializationContext ctxt, TypeDeserializer typeDeserializer) throws IOException, JsonProcessingException {
            return deserialize(jp, ctxt); //To change body of generated methods, choose Tools | Templates.
        }
    }

}
