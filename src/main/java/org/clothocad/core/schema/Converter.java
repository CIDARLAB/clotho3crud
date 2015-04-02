/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import java.util.Map;
import java.util.Set;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
@Slf4j
public abstract class Converter<T extends ObjBase> {

    private Set<Schema> canConvert;
    private Set<String> canConvertNames;
    private Class<T> convertsTo;
    private Schema convertsToSchema;

    public Converter(Schema convertsTo, Set<Schema> canConvert, Set<String> canConvertNames) {
        this.convertsToSchema = convertsTo;
        
        for (String name : canConvertNames) {
            canConvert.add(new InferredSchema(name));
        }
        this.canConvert = canConvert;

        try {
            this.convertsTo = (Class<T>) convertsToSchema.getEnclosedClass(null);
        } catch (ClassNotFoundException ex) {
            log.warn("Couldn't get enclosed class in schema named {}", convertsToSchema.getName());
        }
        this.canConvertNames = canConvertNames;
    }

    public boolean canConvert(Schema type) {
        return (canConvert.contains(type) || canConvertNames.contains(type.getName()));
    }

    public boolean canConvert(String schemaName) {
        return canConvertNames.contains(schemaName);
    }

    public Schema convertsTo() {
        return convertsToSchema;
    }

    public Set<Schema> getCanConvertSchemas() {
        return canConvert;
    }

    public T convert(Map data, Schema type) {
        if (!canConvert(type)) {
            throw new IllegalArgumentException();
        } else {
            try {
                return guardedConvert(data, type);
            } catch (UnsupportedOperationException e) {
                return guardedConvert(data, type.getName());
            }
        }
    }

    public T convert(Map data, String type) {
        if (!canConvert(type)) {
            throw new IllegalArgumentException();
        } else {
            return guardedConvert(data, type);
        }
    }

    protected abstract T guardedConvert(Map data, Schema type);

    protected abstract T guardedConvert(Map data, String schemaName);
}
