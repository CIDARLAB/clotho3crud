/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.annotation.JsonTypeInfo;
import static com.fasterxml.jackson.annotation.JsonTypeInfo.As.EXTERNAL_PROPERTY;
import static com.fasterxml.jackson.annotation.JsonTypeInfo.As.PROPERTY;
import static com.fasterxml.jackson.annotation.JsonTypeInfo.As.WRAPPER_ARRAY;
import static com.fasterxml.jackson.annotation.JsonTypeInfo.As.WRAPPER_OBJECT;
import com.fasterxml.jackson.databind.DeserializationConfig;
import com.fasterxml.jackson.databind.JavaType;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectMapper.DefaultTypeResolverBuilder;
import com.fasterxml.jackson.databind.SerializationConfig;
import com.fasterxml.jackson.databind.cfg.MapperConfig;
import com.fasterxml.jackson.databind.jsontype.NamedType;
import com.fasterxml.jackson.databind.jsontype.TypeDeserializer;
import com.fasterxml.jackson.databind.jsontype.TypeIdResolver;
import com.fasterxml.jackson.databind.jsontype.TypeSerializer;
import com.fasterxml.jackson.databind.jsontype.impl.AsArrayTypeDeserializer;
import com.fasterxml.jackson.databind.jsontype.impl.AsExternalTypeDeserializer;
import com.fasterxml.jackson.databind.jsontype.impl.AsPropertyTypeDeserializer;
import com.fasterxml.jackson.databind.jsontype.impl.AsWrapperTypeDeserializer;
import java.util.Collection;

/**
 *
 * @author spaige
 */
public class WideningDefaultTypeResolverBuilder extends DefaultTypeResolverBuilder {

    public WideningDefaultTypeResolverBuilder(ObjectMapper.DefaultTyping t) {
        super(t);
    }
    public WideningDefaultTypeResolverBuilder() {
        super(ObjectMapper.DefaultTyping.NON_FINAL);
        typeProperty("schema");
        inclusion(JsonTypeInfo.As.PROPERTY);
        
    }

    @Override
    protected TypeIdResolver idResolver(MapperConfig<?> config, JavaType baseType, Collection<NamedType> subtypes, boolean forSer, boolean forDeser) {
        if (_customIdResolver != null) {
            return _customIdResolver;
        }
        if (_idType == null) {
            throw new IllegalStateException("Can not build, 'init()' not yet called");
        }
        if (_idType == JsonTypeInfo.Id.CLASS){
            return new WideningClassNameIdResolver(baseType, config.getTypeFactory());
        }
        return super.idResolver(config, baseType, subtypes, forSer, forDeser); //To change body of generated methods, choose Tools | Templates.
    }

    /*@Override
    public TypeDeserializer buildTypeDeserializer(DeserializationConfig config, JavaType baseType, Collection<NamedType> subtypes) {
        if (!useForType(baseType)) return null;
                if (_idType == JsonTypeInfo.Id.NONE) {
            return null;
        }

        TypeIdResolver idRes = idResolver(config, baseType, subtypes, false, true);
        
        // First, method for converting type info to type id:
        switch (_includeAs) {
        case WRAPPER_ARRAY:
            return new AsArrayTypeDeserializer(baseType, idRes,
                    _typeProperty, _typeIdVisible, _defaultImpl);
        case PROPERTY:
            return new WideningAsPropertyTypeDeserializer(baseType, idRes,
                    _typeProperty, _typeIdVisible, _defaultImpl);
        case WRAPPER_OBJECT:
            return new AsWrapperTypeDeserializer(baseType, idRes,
                    _typeProperty, _typeIdVisible, _defaultImpl);
        case EXTERNAL_PROPERTY:
            return new AsExternalTypeDeserializer(baseType, idRes,
                    _typeProperty, _typeIdVisible, _defaultImpl);
        }
        throw new IllegalStateException("Do not know how to construct standard type serializer for inclusion type: "+_includeAs);
    }*/

    
}
