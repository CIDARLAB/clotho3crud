/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.mongodb.DBCollection;
import com.mongodb.DBObject;
import com.mongodb.ReadPreference;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.jongo.bson.Bson;
import org.jongo.query.Query;
import org.jongo.query.QueryFactory;

/**
 *
 * @author spaige
 */
public class RefFindOne {

    private final ExtendedUnmarshaller unmarshaller;
    private final DBCollection collection;
    private final ReadPreference readPreference;
    private final Query query;
    private Query fields;
    private final QueryFactory queryFactory;

    RefFindOne(DBCollection collection, ReadPreference readPreference, ExtendedUnmarshaller unmarshaller, QueryFactory queryFactory, String query, Object... parameters) {
        this.unmarshaller = unmarshaller;
        this.collection = collection;
        this.readPreference = readPreference;
        this.queryFactory = queryFactory;
        this.query = this.queryFactory.createQuery(query, parameters);
    }

    public <T> T as(final Class<T> clazz) {
        DBObject result = collection.findOne(query.toDBObject(), getFieldsAsDBObject(), readPreference);
        if (result == null) return null;
        
        InstantiatedReferencesQueryingCache cache = new InstantiatedReferencesQueryingCache(queryFactory, unmarshaller, collection);
        
        T resultObj;
        if (ObjBase.class.isAssignableFrom(clazz)) {
            resultObj = (T) cache.makeValue(new ObjectId(result.get("_id")), result, (Class<? extends ObjBase>) clazz);
        } else {
            resultObj  = unmarshaller.unmarshall(Bson.createDocument(result), clazz, cache);
        }

        return resultObj;
    }

    public RefFindOne projection(String fields) {
        this.fields = queryFactory.createQuery(fields);
        return this;
    }

    public RefFindOne projection(String fields, Object... parameters) {
        this.fields = queryFactory.createQuery(fields, parameters);
        return this;
    }

    private DBObject getFieldsAsDBObject() {
        return fields == null ? null : fields.toDBObject();
    }    
}
