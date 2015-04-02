/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
public class Converters {
    
    private Set<Converter> converters = new HashSet<>();
    
    public void addConverter(Converter converter){
        converters.add(converter);
    }

    public Iterable<Schema> getConverterSchemas(Schema schema) {
        //return schemas that can be converted to 'schema' 
        
        Set<Schema> results = new HashSet<>();
        for (Converter converter : converters){
            if (converter.convertsTo().equals(schema)) results.addAll(converter.getCanConvertSchemas());
        }
        
        return results;
    }

    public Converter getConverter(Schema from, Schema to){
        //TODO: graph search for multi-step conversion
        //TODO: CompositeConvereter
        
        
        for (Converter converter : converters){
            if (converter.convertsTo().equals(to) && converter.canConvert(from)){
                return converter;
            }
        }
        
        return null;
    }
}
