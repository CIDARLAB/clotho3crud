/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.JsonToken;
import com.fasterxml.jackson.core.util.JsonParserSequence;
import com.fasterxml.jackson.databind.BeanProperty;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.JavaType;
import com.fasterxml.jackson.databind.JsonDeserializer;
import com.fasterxml.jackson.databind.jsontype.TypeIdResolver;
import com.fasterxml.jackson.databind.jsontype.impl.AsPropertyTypeDeserializer;
import com.fasterxml.jackson.databind.util.TokenBuffer;
import java.io.IOException;

/**
 *
 * @author spaige
 */
public class WideningAsPropertyTypeDeserializer extends AsPropertyTypeDeserializer {

    public WideningAsPropertyTypeDeserializer(AsPropertyTypeDeserializer src, BeanProperty property) {
        super(src, property);
    }

    public WideningAsPropertyTypeDeserializer(JavaType bt, TypeIdResolver idRes, String typePropertyName, boolean typeIdVisible, Class<?> defaultImpl) {
        super(bt, idRes, typePropertyName, typeIdVisible, defaultImpl);
    }

    @Override
    public Object deserializeTypedFromObject(JsonParser jp, DeserializationContext ctxt) throws IOException, JsonProcessingException {
        // but first, sanity check to ensure we have START_OBJECT or FIELD_NAME
        JsonToken t = jp.getCurrentToken();
        if (t == JsonToken.START_OBJECT) {
            t = jp.nextToken();
        } else if (t == JsonToken.START_ARRAY) {
            /* This is most likely due to the fact that not all Java types are
             * serialized as JSON Objects; so if "as-property" inclusion is requested,
             * serialization of things like Lists must be instead handled as if
             * "as-wrapper-array" was requested.
             * But this can also be due to some custom handling: so, if "defaultImpl"
             * is defined, it will be asked to handle this case.
             */
            return _deserializeTypedUsingDefaultImpl(jp, ctxt, null);
        } else if (t != JsonToken.FIELD_NAME) {
            return _deserializeTypedUsingDefaultImpl(jp, ctxt, null);
        }
        // Ok, let's try to find the property. But first, need token buffer...
        TokenBuffer tb = null;

        for (; t == JsonToken.FIELD_NAME; t = jp.nextToken()) {
            String name = jp.getCurrentName();
            jp.nextToken(); // to point to the value
            if (_typePropertyName.equals(name)) { // gotcha!
                return modifiedDeserializeTypedForId(jp, ctxt, tb);
            }
            if (tb == null) {
                tb = new TokenBuffer(null);
            }
            tb.writeFieldName(name);
            tb.copyCurrentStructure(jp);
        }
        return _deserializeTypedUsingDefaultImpl(jp, ctxt, tb);
    }

   protected Object modifiedDeserializeTypedForId(JsonParser jp, DeserializationContext ctxt,
            TokenBuffer tb)
        throws IOException, JsonProcessingException
    {
        String typeId = jp.getText();
        JsonDeserializer<Object> deser = modifiedFindDeserializer(ctxt, typeId);
        if (_typeIdVisible) { // need to merge id back in JSON input?
            if (tb == null) {
                tb = new TokenBuffer(null);
            }
            tb.writeFieldName(jp.getCurrentName());
            tb.writeString(typeId);
        }
        if (tb != null) { // need to put back skipped properties?
            jp = JsonParserSequence.createFlattened(tb.asParser(jp), jp);
        }
        // Must point to the next value; tb had no current, jp pointed to VALUE_STRING:
        jp.nextToken(); // to skip past String value
        // deserializer should take care of closing END_OBJECT as well
        return deser.deserialize(jp, ctxt);
    }
    
    protected JsonDeserializer<Object> modifiedFindDeserializer(DeserializationContext ctxt,
            String typeId)
        throws IOException, JsonProcessingException
    {
        JsonDeserializer<Object> deser;

        synchronized (_deserializers) {
            deser = _deserializers.get(typeId);
            if (deser == null) {
                JavaType type = _idResolver.typeFromId(typeId);
                if (type == null) {
                    // As per [JACKSON-614], use the default impl if no type id available:
                    if (_defaultImpl == null) {
                        throw ctxt.unknownTypeException(_baseType, typeId);
                    }
                    deser = _findDefaultImplDeserializer(ctxt);
                } else {
                    /* 16-Dec-2010, tatu: Since nominal type we get here has no (generic) type parameters,
                     *   we actually now need to explicitly narrow from base type (which may have parameterization)
                     *   using raw type.
                     *   
                     *   One complication, though; can not change 'type class' (simple type to container); otherwise
                     *   we may try to narrow a SimpleType (Object.class) into MapType (Map.class), losing actual
                     *   type in process (getting SimpleType of Map.class which will not work as expected)
                     */
                    if (_baseType != null && _baseType.getClass() == type.getClass()) {
                        try{ //the actual modification behind layers of final methods!
                            type = _baseType.narrowBy(type.getRawClass());
                        } catch (IllegalArgumentException iae){
                            type = _baseType.widenBy(type.getRawClass());
                        }
                    }
                    deser = ctxt.findContextualValueDeserializer(type, _property);
                }
                _deserializers.put(typeId, deser);
            }
        }
        return deser;
    }
    
}
