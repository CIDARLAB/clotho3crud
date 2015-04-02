 /*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.persistence;



import groovy.lang.Tuple;
import java.net.UnknownHostException;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;


/**
 *
 * @author J. Christopher Anderson
 * @author Douglas Densmore
 */
public interface ClothoConnection {

    /**
     * Returns true if the database backing this connection adheres to the
     * Clotho data model.
     * XXX: what does it return if there is no connection?
     * @return
     */
    boolean isAClothoDatabase();

    /**
     * Programmatically tell the connection to disconnect itself
     * from its database
     * 
     * throws XXX if errors
     */
    void disconnect();

    /**
     * Determine whether the database is connected currently
     * @return true if currently connected
     */
    boolean isConnected();

    /**
     * Saves the given object and its references to the database.
     * 
     * throws XXX if errors
     * 
     * @param obj
     */
    void save(ObjBase obj);
    void save(Map obj);
    
    /**
     * Saves the given collection of objects to the database.
     * @param objs
     * @return the number of objects successfully saved
     */
    void saveAll(Iterable<ObjBase> objs);
    int saveBSON(Collection<Map> objs);

    
    /**
     * Delete the object from the database.
     * @param obj
     * @return
     */
    void delete( ObjBase obj );
    void delete(ObjectId id);

    /**
     * Deletes the given set of objects from the database.
     * @param objs
     * @return number of objects deleted
     */
    int delete( Collection<ObjBase> objs );
    

    /**
     * Returns the time the given ObjBase object was modified in the database.
     * @param obj
     * @return
     */
    Date getTimeModified( ObjBase obj );

    /**
     * Gets the object with the given uuid as the specified class,
     * or null if the object does not exist in the database.
     * 
     * @param type
     * @param uuid
     * @return
     */
    <T extends ObjBase> T get(Class<T> type, ObjectId uuid);
    Map<String,Object> getAsBSON(ObjectId uuid);
    Map<String,Object> getAsBSON(ObjectId uuid, Set<String> fields);
    

    /**
     * Search the database using a query using MongoDB semantics.
     * See: 
     * http://docs.mongodb.org/manual/core/read-operations/
     * http://docs.mongodb.org/manual/reference/operator/
     * 
     * Objects will be deserialized to their 'default' type (the type they were first saved as.)
     * 
     * Type information is stored in the 'type' field as the canonical java name of the default type
     * @param query
     * @return
     */
    List<ObjBase> get(Map query);
    List<ObjBase> get(Map query, int hitmax);
    
    <T extends ObjBase> List<T> get(Class<T> type, Map query);
    <T extends ObjBase> List<T> get(Class<T> type, Map query, int hitmax);
    
    List<Map<String,Object>> getAsBSON(Map query);
    List<Map<String,Object>> getAsBSON(Map query, int hitmax, Set<String> fields);

    <T extends ObjBase> T getOne(Class<T> type, Map query);

    Map<String,Object> getOneAsBSON(Map query);
    Map<String,Object> getOneAsBSON(Map query, Set<String> fields);
    
    <T extends ObjBase> List<T> getAll(Class<T> type);
    
    //Deletes everything
    void deleteAll();

    public boolean exists(ObjectId id);
    
    public List<Map> getCompletionData();
 
}
