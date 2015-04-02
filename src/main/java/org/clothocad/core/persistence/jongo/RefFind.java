/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.ReadPreference;
import com.thoughtworks.proxy.toys.hotswap.Swappable;
import java.util.ArrayList;
import java.util.List;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.ObjBase;
import org.jongo.bson.Bson;
import org.jongo.query.Query;
import org.jongo.query.QueryFactory;

/**
 * Copy of Jongo Find, for modification purposes
 *
 * @author spaige
 */
@Slf4j
public class RefFind {

    private final DBCollection collection;
    private final ReadPreference readPreference;
    private final ExtendedUnmarshaller unmarshaller;
    private final QueryFactory queryFactory;
    private final Query query;
    private Query fields, sort, hint;
    private Integer limit, skip;

    /**
     *
     * @param collection
     * @param readPreference
     * @param unmarshaller
     * @param queryFactory
     * @param query
     * @param parameters
     */
    public RefFind(DBCollection collection, ReadPreference readPreference, ExtendedUnmarshaller unmarshaller, QueryFactory queryFactory, String query, Object... parameters) {
        this.readPreference = readPreference;
        this.unmarshaller = unmarshaller;
        this.collection = collection;
        this.queryFactory = queryFactory;
        this.query = this.queryFactory.createQuery(query, parameters);
    }

    //Greedy instead of lazy resolution of results - compare to jongo Find
    public <T> Iterable<T> as(final Class<T> clazz) {
        DBCursor cursor = new DBCursor(collection, query.toDBObject(), getFieldsAsDBObject(), readPreference);
        addOptionsOn(cursor);
        List<T> out = new ArrayList<>();
        InstantiatedReferencesCache cache = new InstantiatedReferencesCache();
        if (ObjBase.class.isAssignableFrom(clazz)) {
            while (cursor.hasNext()) {
                DBObject next = cursor.next();
                String id = unmarshaller.getId(Bson.createDocument(next));
                ObjBase current;
                if (cache.has(id)) {
                    //get pre-existing value, update it
                    current = (ObjBase) unmarshaller.unmarshall(Bson.createDocument(next), cache.get(id), cache);
                } else {
                    //not in cache, deserialize and add to cache
                    current = (ObjBase) unmarshaller.unmarshall(Bson.createDocument(next), clazz, cache);
                    //TODO: remove 'done' values from cache processing queue
                    cache.addValue(current.getId().toString(), current);
                }
                out.add((T) current);
            }
            // resolve any referenced objects
            while (!cache.done()) {
                ObjBase current = cache.getNextUndone();
                Query currentQuery = queryFactory.createQuery("{_id:#}", current.getId().toString());
                DBObject result = collection.findOne(currentQuery.toDBObject(), null, readPreference);
                if (result == null) {
                    log.warn("couldn't find referenced object '{}'", current.getId());
                    continue;
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

        return out;
    }

    private void addOptionsOn(DBCursor cursor) {
        if (limit != null) {
            cursor.limit(limit);
        }
        if (skip != null) {
            cursor.skip(skip);
        }
        if (sort != null) {
            cursor.sort(sort.toDBObject());
        }
        if (hint != null) {
            cursor.hint(hint.toDBObject());
        }
    }

    public RefFind projection(String fields) {
        this.fields = queryFactory.createQuery(fields);
        return this;
    }

    public RefFind projection(String fields, Object... parameters) {
        this.fields = queryFactory.createQuery(fields, parameters);
        return this;
    }

    public RefFind limit(int limit) {
        this.limit = limit;
        return this;
    }

    public RefFind skip(int skip) {
        this.skip = skip;
        return this;
    }

    public RefFind sort(String sort) {
        this.sort = queryFactory.createQuery(sort);
        return this;
    }

    public RefFind hint(final String hint) {
        this.hint = queryFactory.createQuery(hint);
        return this;
    }

    private DBObject getFieldsAsDBObject() {
        return fields == null ? null : fields.toDBObject();
    }
}
