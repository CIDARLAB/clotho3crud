/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import org.clothocad.core.persistence.jackson.JSONFilter;
import static com.fasterxml.jackson.annotation.JsonInclude.Include.NON_NULL;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.type.TypeReference;
import static com.fasterxml.jackson.databind.DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES;
import com.fasterxml.jackson.databind.ObjectMapper;
import static com.fasterxml.jackson.databind.SerializationFeature.FAIL_ON_EMPTY_BEANS;
import com.mongodb.BasicDBObject;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.MongoClient;
import com.mongodb.WriteConcernException;
import static com.mongodb.MongoException.DuplicateKey;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.inject.Inject;
import javax.inject.Named;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.authz.Permission;
import org.apache.shiro.authz.permission.RolePermissionResolver;
import org.apache.shiro.authz.permission.WildcardPermission;
import org.bson.BSONObject;
import org.bson.LazyBSONList;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.AuthGroup;
import org.clothocad.core.security.ClothoAccount;
import org.clothocad.core.security.ClothoAction;
import org.clothocad.core.security.CredentialStore;
import org.clothocad.core.security.PermissionsOnObject;
import org.clothocad.core.util.JSON;
import org.jongo.MongoCollection;
import org.jongo.ResultHandler;
import org.python.google.common.collect.Lists;

/**
 *
 * @author spaige
 */
@Slf4j
public class JongoConnection implements ClothoConnection, CredentialStore, RolePermissionResolver {

    protected RefJongo jongo;
    protected DB db;
    protected MongoClient client;
    protected ObjectMapper mapper;
    //TODO: split into separate collections per top-level schema
    protected RefMongoCollection data;
    protected MongoCollection cred;
    protected MongoCollection roles;
    protected DBCollection rawDataCollection;
    
    protected DBClassLoader classLoader;
    private static final TypeReference<Map<String, Object>> STRINGMAP = new TypeReference<Map<String, Object>>() {
    };

    @Inject
    public JongoConnection(@Named("dbport") int port, @Named("dbhost") String host, @Named("dbname") String dbname, DBClassLoader dbClassLoader) throws UnknownHostException {
        log.info("Using mongo database '{}@{}:{}'", dbname, host, port);
        
        db = new MongoClient(host, port).getDB(dbname);
        rawDataCollection = db.getCollection("data");
        classLoader = dbClassLoader;
        connect();
    }

    private void connect() throws UnknownHostException {
        //TODO: cover reconnect case?

        //Mimic Jongo customization         
        mapper = new ObjectMapper();
        mapper.disable(FAIL_ON_UNKNOWN_PROPERTIES);
        mapper.setSerializationInclusion(NON_NULL);
        mapper.disable(FAIL_ON_EMPTY_BEANS);
        //mapper.enableDefaultTyping(ObjectMapper.DefaultTyping.OBJECT_AND_NON_CONCRETE);
        //jongo mimicking over
        mapper.registerModule(new JSON.ClothoJacksonModule());

        /**
         * redundant with ObjBase annotations
         *
         * mapper.disable(AUTO_DETECT_GETTERS);
         * mapper.disable(AUTO_DETECT_IS_GETTERS);
         **/

        jongo = new RefJongo(db, new ClothoMapper());
        data = jongo.getCollection("data");
        cred = jongo.getCollection("cred");
        roles = jongo.getCollection("roles");
    }

    //Do we really need this?
    @Override
    public boolean isAClothoDatabase() {
        //TODO
        return db != null;
    }

    @Override
    public void disconnect() {
        client.close();
    }

    @Override
    public boolean isConnected() {
        return db != null;
    }

    @Override
    public void save(ObjBase obj) {
        try{
            data.save(obj);
        }
        catch (DuplicateKey e){
            data.update("{_id:#}", obj.getId().toString()).with(obj);
        }
    }

    @Override
    public void save(Map obj) {
        obj = mongifyIdField(obj);
        if (obj.get("_id") == null){
            //must create new id
            rawDataCollection.insert(new BasicDBObject(obj));
            return;
        }
        DBObject idQuery = new BasicDBObject("_id", obj.get("_id"));
        obj.remove("_id");
        Map<String,Object> setExpression = new HashMap<>();
        setExpression.put("$set", obj);
        try {
            rawDataCollection.update(idQuery, new BasicDBObject(setExpression),true, false);  //upsert true, multi false        
        } catch (WriteConcernException ex) {
            log.error("Invalid JSON/failed to save object", obj.toString());
        }
    }

    @Override
    public void saveAll(Iterable<ObjBase> objs) {
        for (ObjBase o : objs) {
            save(o);
        }
    }

    @Override
    public int saveBSON(Collection<Map> objs) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void delete(ObjBase obj) {
        //TODO: check for references in database
        //TODO: move to 'deleted' collection instead of hard deleting
        delete(obj.getId());
    }

    @Override
    public void delete(ObjectId id) {
        data.remove("{_id:#}", id.toString());
    }

    @Override
    public int delete(Collection<ObjBase> objs) {
        int i = 0;
        for (ObjBase o : objs) {
            try {
                delete(o);
                i++;
            } catch (Exception e) {
                log.error("Error while deleting object in collection", e);
                //aggregate errors and pass back to user - not sure on details yet 
            }
        }
        return i;
    }

    @Override
    public Date getTimeModified(ObjBase obj) {
        //TODO: just fetch LastModified field instead of /entire object
        // use data.projection('{lastModified:1}').as(...?
        ObjBase result = data.resolvingFindOne("{_id:#}", obj.getId().toString()).as(ObjBase.class);
        return result.getLastModified();
    }

    @Override
    public <T extends ObjBase> T get(Class<T> type, ObjectId uuid) {
        bindClassLoader();
        return data.resolvingFindOne("{_id:#}", uuid.toString()).as(type);
    }

    @Override
    public Map<String, Object> getAsBSON(ObjectId uuid) {
        return getAsBSON(uuid, null);
    }
    
    @Override
    public Map<String, Object> getAsBSON(ObjectId uuid, Set<String> filter) {
        return data.findOne("{_id:#}", uuid.toString()).projection(generateProjection(filter)).map(DemongifyHandler.get());
    }

    @Override
    public List<ObjBase> get(Map query) {
        return get(query, Persistor.SEARCH_MAX);
    }

    
    @Override
    public List<ObjBase> get(Map query, int hitmax) {
        bindClassLoader();
        return Lists.newArrayList(data.resolvingFind(serialize(mongifyIdField(query))).limit(hitmax).as(ObjBase.class));
    }

    @Override
    public <T extends ObjBase> List<T> get(Class<T> type, Map query) {
        return get(type, query, Persistor.SEARCH_MAX);
    }
    
    
    @Override
    public <T extends ObjBase> List<T> get(Class<T> type, Map query, int hitmax) {
        bindClassLoader();
        return Lists.newArrayList(data.resolvingFind(serialize(mongifyIdField(query))).limit(hitmax).as(type));
    }

    @Override
    public List<Map<String, Object>> getAsBSON(Map query) {    
        return getAsBSON(query, Persistor.SEARCH_MAX, null);
    }
    
    
    @Override
    public List<Map<String, Object>> getAsBSON(Map query, int hitmax, Set<String> filter) {
        return Lists.newArrayList(data.find(serialize(mongifyIdField(query)))
                .limit(hitmax).projection(generateProjection(filter))
                .map(DemongifyHandler.get()));
    }

    @Override
    public <T extends ObjBase> T getOne(Class<T> type, Map query) {
        bindClassLoader();
        return data.resolvingFindOne(serialize(mongifyIdField(query))).as(type);
    }

    @Override
    public Map<String, Object> getOneAsBSON(Map query) {
        return getOneAsBSON(query, null);
    }
    @Override
    public Map<String, Object> getOneAsBSON(Map query, Set<String> filter) {
        return data.findOne(serialize(mongifyIdField(query))).projection(generateProjection(filter)).map(DemongifyHandler.get());
    }

    @Override
    public <T extends ObjBase> List<T> getAll(Class<T> type) {
        return Lists.newArrayList(data.resolvingFind("{schema:#}", type.getCanonicalName()).as(type));
    }

    @Override
    public void deleteAll() {
        rawDataCollection.drop();
    }

    @Override
    public boolean exists(ObjectId id) {
        return rawDataCollection.find(new BasicDBObject("_id", id.toString())).hasNext();
    }

    @Override
    public ClothoAccount getAccount(String username) {
        return cred.findOne("{_id:#}", username).as(ClothoAccount.class);
        }

    @Override
    public void saveAccount(ClothoAccount account) {
        cred.save(account);
    }

    public AuthGroup getGroup(String groupName){
        return roles.findOne("{_id:#}", groupName).as(AuthGroup.class);
    }

    public void saveGroup (AuthGroup group){
        roles.save(group);
    }

    @Override
    public void deleteAllCredentials() {
        cred.drop();
        roles.drop();
    }
    
    @Override
    public List<Map> getCompletionData(){
        DBCursor cursor = rawDataCollection.find();
        Iterator<DBObject> iter = cursor.iterator();
        List<Map> out = new ArrayList<Map>();
        int i = 0;
        while(iter.hasNext()){
            try {
                DBObject temp = iter.next();
                Map map = new HashMap();
                map.put("name", temp.get("name"));
                map.put("schema", temp.get("schema"));
                map.put("id", temp.get("_id").toString());
                if(temp.containsKey("description")) {
                    map.put("description", temp.get("description"));
                } else if(temp.containsKey("shortDescription")) {
                    map.put("description", temp.get("shortDescription"));
                }
                out.add(map);
            } catch(Exception err) {
                err.printStackTrace();
            }
            i++;
        }
        return out;
    }
    private List<Map<String, Object>> mappify(Iterable<JSONFilter> objs) {
        List<Map<String, Object>> out = new ArrayList<>();
        for (JSONFilter obj : objs) {
            out.add(mappify(obj));
        }
        return out;
    }
    
    private Map<String, Object> mappify(JSONFilter obj) {

        return mapper.convertValue(obj, STRINGMAP);
    }

    private String serialize(Object obj) {
        try {
            return mapper.writeValueAsString(obj);
        } catch (JsonProcessingException ex) {
            throw new RuntimeException(ex);
        }
    }

    //renames "id" key to "_id", and replaces ObjectId with String
    protected static Map<String,Object> mongifyIdField(Map<String, Object> obj) {
        if (obj.containsKey(ID)){
            Map<String,Object> copy = new HashMap<>();
            copy.putAll(obj);
            obj = copy;
            Object id = obj.get(ID);
            if (id instanceof ObjectId){
                id = ((ObjectId) id).toString();
            }
            obj.remove(ID);
            obj.put("_id", id);
        }
        return obj;
    }
    
    //renames "_id" to "id" 
    //doesn't recurse on sub-objects, since no embedded object will have an objectid
    //mutates instead of returning new object
    protected static Map<String,Object> demongifyIdField(Map<String,Object> obj){
        if (obj.containsKey("_id")){
            Object id = obj.get("_id");
            obj.remove("_id");
            obj.put(ID, id);
        }
        return obj;
    }

    private void bindClassLoader() {
        //Set up classloader to look in db for class definitions.
        // this may not be the right place to do this - not sure if better to do this on thread spawn always
        //classloader is singleton that delegates to root classloader, so no classloader hairiness should happen
        Thread.currentThread().setContextClassLoader(classLoader);
    }
    
    private String generateProjection(Set<String> fields){
        if (fields == null) return "{}";
        StringBuilder builder = new StringBuilder();
        builder.append("{");
        for (String field: fields){
            builder.append(String.format("%s:1, ", field));
        }
        builder.append("}");
        
        return builder.toString();
    }
    
    @Override
    public Collection<Permission> resolvePermissionsInRole(String role) {
        Collection<Permission> permissions = new HashSet<>();
        AuthGroup group = roles.findOne("{_id:#}", role).as(AuthGroup.class);
        
        if (group == null) return permissions;
        
        for (PermissionsOnObject p : group.getPermissions().values()){
            for (String permString : p.getPermissions()){
                permissions.add(new WildcardPermission(permString));
            }
        }
        return permissions;
    }
    
    @Override
    public Map<String, Set<ClothoAction>> getUserPermissions(ObjectId id){
        Map<String, Set<ClothoAction>> permissionsByUser = new HashMap<>();
        
        for (ClothoAccount account : cred.find("{'authzInfo.permissions.#':{$exists:true}}", id)
            .projection("{'authzInfo.permissions.#':1, '@class':1}", id).as(ClothoAccount.class)){
            permissionsByUser.put(account.getId(), account.getActions(id));
        }
        
        return permissionsByUser;
    }

    @Override
    public Map<String, Set<ClothoAction>> getGroupPermissions(ObjectId id){
        Map<String, Set<ClothoAction>> permissionsByGroup = new HashMap<>();
        
        for (AuthGroup group : roles.find("{'permissions.#':{$exists:true}}", id)
            .projection("{'permissions.#':1, '$class':1}", id).as(AuthGroup.class)){
            permissionsByGroup.put(group.getName(), group.getActions(id));
        }
        
        return permissionsByGroup;        
    }
    
    protected static class DemongifyHandler implements ResultHandler<Map<String,Object>>{
        //this will have problems if there are circular dstructures, though that should be impossible
        @Override
        public Map<String,Object> map(DBObject result) {
            Map<String,Object> resultMap = toMap(result);
            //recurse on any sub objects
            for (String key : resultMap.keySet()){
                 Object value = resultMap.get(key);
                 if (value instanceof DBObject){
                     resultMap.put(key,toMapOrList((DBObject)value));
                 }
            }
            return demongifyIdField(resultMap);
        }
        
        //some dbobjects don't support toMap
        private Map<String,Object> toMap(BSONObject dbObject){
            BasicDBObject basicResult = new BasicDBObject();
            basicResult.putAll(dbObject);
            Map<String,Object> resultMap = basicResult.toMap();
            return resultMap;
        }
        
        private Object toMapOrList(DBObject value) {
            if (value instanceof LazyBSONList) {
                //convert members
                List convertedList = new ArrayList();
                for (Object element : (List) value) {
                    if (element instanceof DBObject) {
                        convertedList.add(toMapOrList((DBObject) element));
                    } else {
                        convertedList.add(element);
                    }
                }
                return convertedList;

            } else {
                return map(value);
            }
        }

        private static DemongifyHandler instance;
        
          public static DemongifyHandler get(){
            if (instance == null){
                instance = new DemongifyHandler();
            }
            return instance;
        }

    } 
}
