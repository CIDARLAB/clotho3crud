/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.mongodb.DB;
import com.mongodb.DBCollection;
import org.jongo.Jongo;
import org.jongo.MongoCollection;
import org.jongo.bson.BsonDBDecoder;
import org.jongo.bson.BsonDBEncoder;

/**
 * Subclassing Jongo so we have hooks to put our own versions of Find, FindOne in
 * 
 * @author spaige
 */
public class RefJongo extends Jongo {

    DB database;
    ClothoMapper mapper;
    
    public RefJongo(DB database) {
        super(database);
        this.database = database;
    }

    public RefJongo(DB database, ClothoMapper mapper) {
        super(database, mapper);
        this.database = database;
        this.mapper = mapper;
    }

    @Override
    public RefMongoCollection getCollection(String name) {
        DBCollection dbCollection = database.getCollection(name);
        dbCollection.setDBDecoderFactory(BsonDBDecoder.FACTORY);
        dbCollection.setDBEncoderFactory(BsonDBEncoder.FACTORY);
        return new RefMongoCollection(dbCollection, mapper);
    }
    
    
}
