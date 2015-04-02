/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.mongodb.DBCollection;
import org.jongo.MongoCollection;

/**
 *
 * @author spaige
 */
public class RefMongoCollection extends MongoCollection {

    ClothoMapper mapper;
    
    
    public RefMongoCollection(DBCollection dbCollection, ClothoMapper mapper) {
        super(dbCollection, mapper);
        this.mapper = mapper;
    }

    public RefFind resolvingFind(String query, Object... parameters){
        DBCollection coll= getDBCollection();
        return new RefFind(coll, coll.getReadPreference(), mapper.getUnmarshaller(), mapper.getQueryFactory(), query, parameters);
    }
    
    public RefFindOne resolvingFindOne(String query, Object... parameters){
        DBCollection coll = getDBCollection();
        return new RefFindOne(coll, coll.getReadPreference(), mapper.getUnmarshaller(), mapper.getQueryFactory(), query, parameters);
    }
    
}
