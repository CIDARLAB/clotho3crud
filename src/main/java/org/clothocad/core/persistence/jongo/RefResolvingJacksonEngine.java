/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.InjectableValues;
import java.io.IOException;
import org.clothocad.core.persistence.jackson.JSONFilter;
import org.clothocad.core.persistence.jackson.JSONViews;
import org.jongo.bson.BsonDocument;
import org.jongo.marshall.MarshallingException;
import org.jongo.marshall.jackson.JacksonEngine;
import org.jongo.marshall.jackson.configuration.Mapping;

/**
 *
 * @author spaige
 */
public class RefResolvingJacksonEngine extends JacksonEngine implements ExtendedUnmarshaller{

    private Mapping mapping;
    public RefResolvingJacksonEngine(Mapping mapping) {
        super(mapping);
        this.mapping = mapping;
    }

    //unmarshalls document to an instance of clazz, using a cache of injectable values
    @Override
    public <T> T unmarshall(BsonDocument document, Class<T> clazz, InjectableValues cache) {
        try {
            return mapping.getReader(clazz).with(cache).readValue(document.toByteArray(), 0, document.getSize());
        } catch (IOException e) {
            String message = String.format("Unable to unmarshall result to %s from content %s", clazz, document.toString());
            throw new MarshallingException(message, e);
        }
    }

    //updates toUpdate using the data in document and a cache of injectable values
    @Override
    public <T> T unmarshall(BsonDocument document, T toUpdate, InjectableValues cache) {
        try{
            return mapping.getReader(toUpdate.getClass()).with(cache).withValueToUpdate(toUpdate).readValue(document.toByteArray(), 0, document.getSize());
        } catch (IOException e) {
            String message = String.format("Unable to unmarshall result to %s from content %s", toUpdate.getClass(), document.toString());
            throw new MarshallingException(message, e);
        }
    }

    @Override
    public String getId(BsonDocument document) {
        try {
            JSONFilter filter = mapping.getReader(JSONFilter.class).withView(JSONViews.IdOnly.class).readValue(document.toByteArray(), 0, document.getSize());
            return filter.getId().toString();
        } catch (IOException e) {
            String message = String.format("Unable to extract id from content %s", document.toString());
            throw new MarshallingException(message, e);
        }
    }
}
