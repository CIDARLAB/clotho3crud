/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.mongodb.DBCollection;
import com.mongodb.DBObject;
import com.mongodb.ReadPreference;
import com.thoughtworks.proxy.toys.hotswap.Swappable;
import org.clothocad.core.datums.ObjBase;
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
        
        InstantiatedReferencesCache cache = new InstantiatedReferencesCache();
        
        T resultObj = unmarshaller.unmarshall(Bson.createDocument(result), clazz, cache);
        if (ObjBase.class.isAssignableFrom(clazz)) {
            ObjBase current = (ObjBase) resultObj;
            cache.addValue(current.getId().toString(), current);
            while (!cache.done()) {
                current = cache.getNextUndone();
                Query currentQuery = queryFactory.createQuery("{_id:#}", current.getId().toString());
                result = collection.findOne(currentQuery.toDBObject(), null, readPreference);
                if (result == null){
                    throw new RuntimeException("Unresolvable reference: No object with id "+ current.getId() +" found.");
                }
                if (current instanceof Swappable){
                    //need to resolve proxy to actual value
                    ObjBase actual = unmarshaller.unmarshall(Bson.createDocument(result), ObjBase.class , cache);
                    ((Swappable) current).hotswap(actual);
                }
                else {
                    //update instance with results of query
                    unmarshaller.unmarshall(Bson.createDocument(result), current, cache);
                }
            }
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
