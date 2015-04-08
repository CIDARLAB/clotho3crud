/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.ReadPreference;
import java.util.ArrayList;
import java.util.List;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
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

    public <T> Iterable<T> as(final Class<T> clazz) {
        DBCursor cursor = new DBCursor(collection, query.toDBObject(), getFieldsAsDBObject(), readPreference);
        addOptionsOn(cursor);
        List<T> out = new ArrayList<>();
        InstantiatedReferencesQueryingCache cache = new InstantiatedReferencesQueryingCache(queryFactory, unmarshaller, collection);
        while (cursor.hasNext()) {
            DBObject next = cursor.next();
            T current;
            if (ObjBase.class.isAssignableFrom(clazz)) {
                current = (T) cache.makeValue(new ObjectId(next.get("_id")), next, (Class<? extends ObjBase>) clazz);
            } else {
                current = unmarshaller.unmarshall(Bson.createDocument(next), clazz, cache);
            }
            out.add((T) current);
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
