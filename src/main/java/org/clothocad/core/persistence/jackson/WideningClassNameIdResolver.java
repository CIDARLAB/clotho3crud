/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.databind.JavaType;
import com.fasterxml.jackson.databind.jsontype.impl.ClassNameIdResolver;
import com.fasterxml.jackson.databind.type.TypeFactory;
import com.fasterxml.jackson.databind.util.ClassUtil;

/**
 *
 * @author spaige
 */
public class WideningClassNameIdResolver extends ClassNameIdResolver {

    public WideningClassNameIdResolver(JavaType baseType, TypeFactory typeFactory) {
        super(baseType, typeFactory);
    }

    @Override
    public JavaType typeFromId(String id) {
        try{
            return super.typeFromId(id);
        } catch (IllegalArgumentException iae){
            try {
                //try the widening type construction instead
                Class<?> cls =  ClassUtil.findClass(id);
                if (cls.isAssignableFrom(_baseType.getRawClass()) || _baseType.getRawClass().equals(JSONFilter.class)) return _baseType;
                //TODO: make this message not suck
                else throw new IllegalArgumentException("Invalid type id '"+id+"' (for id type 'Id.class'): ");
            } catch (ClassNotFoundException cnfe){
                throw new IllegalArgumentException("Invalid type id '"+id+"' (for id type 'Id.class'): no such class found");
            } 
            catch (Exception e){
                //could neither narrow or widen that type combination
                throw new IllegalArgumentException("Invalid type id '"+id+"' (for id type 'Id.class'): "+e.getMessage(), e);
            } 
            
        }
    }
    
    
}
