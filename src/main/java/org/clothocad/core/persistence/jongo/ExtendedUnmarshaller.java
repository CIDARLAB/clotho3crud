/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.InjectableValues;
import com.mongodb.DBObject;
import org.jongo.bson.BsonDocument;
import org.jongo.marshall.Unmarshaller;

/**
 *
 * @author spaige
 */
public interface ExtendedUnmarshaller extends Unmarshaller {
    <T> T unmarshall(BsonDocument document, Class clazz, InjectableValues cache);
    <T> T unmarshall(BsonDocument document, T toUpdate, InjectableValues cache);

    public String getId(BsonDocument document);
}
