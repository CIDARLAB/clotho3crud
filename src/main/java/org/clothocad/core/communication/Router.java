package org.clothocad.core.communication;

import com.fasterxml.jackson.core.type.TypeReference;
import com.google.common.collect.Lists;
import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.inject.Inject;
import javax.inject.Singleton;
import lombok.Getter;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.Subject;
import static org.clothocad.core.communication.Channel.create;
import static org.clothocad.core.communication.Channel.destroy;
import static org.clothocad.core.communication.Channel.log;
import static org.clothocad.core.communication.Channel.submit;

import org.clothocad.core.datums.ObjectId;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.util.JSON;

@Slf4j
@Singleton
public class Router {

    @Getter
    protected Persistor persistor;
    
    @Getter
    protected ClothoRealm realm;
    
    @Inject
    public Router(Persistor persistor, ClothoRealm realm) {
        minds = new HashMap<>();
        this.persistor = persistor;
        this.realm = realm;
    }

    // send message    
    public void sendMessage(ClientConnection connection, Message message) {
        try {
            log.debug(JSON.serialize(message));
        } catch (IOException e) {
            log.debug("failed to serialize message: {}", message);
        }
        connection.send(message);
    }

    // receive message
    public void receiveMessage(ClientConnection connection, Message request) {

        //bind context to request
        Subject subject = SecurityUtils.getSubject();
        Mind mind;
        boolean wasAuthenticated = false;
        if (subject.isAuthenticated()){
            wasAuthenticated = true;
            mind = getAuthenticatedMind(subject.getPrincipal().toString(), connection, request);
            mind.setConnection(connection);
        } else {
            mind = getMind(connection);
            //Find a better way to become anonymous user?
            subject.login(ClothoRealm.getAnonymousUserToken());
        }
        ServerSideAPI api = new ServerSideAPI(mind, persistor, this, realm, request.getRequestId(), new MessageOptions(request.getOptions()));

        Object data = request.getData();
        
        Object response = null;
        try {
            switch (request.getChannel()) {                
                case autocomplete:
                    response = api.autocomplete(data.toString());
                    break;
                case submit:
                    response = api.submit(data);
                    break;
                case clear:
                    api.clear();
                    break;
                case login:
                    Map map = (Map) data;
                    if (!wasAuthenticated){
                        //log out of anonymous user
                        //XXX: all this is madness
                        subject.logout();
                    }
                    response = api.login(map.get("username").toString(), map.get("credentials").toString());
                    //we only reach this point if login was successful
                    if (!wasAuthenticated){
                        //remove the 'anonymous' mind
                        //currently this means you lose environment state if you login
                        //we could do something more sophisticated like merge the anonymous environment and the persisted mind, but that could get complicated
                        minds.remove(connection.getId());
                    }
                    break;
                case updatePassword:
                    Map updatedPassword = (Map) data;
                    response  = api.changePassword(updatedPassword.get("password").toString());
                    
                    break;
                
                case createUser:
                    Map newusermap = (Map) data;
                    
                    if(newusermap.containsKey("password"))
                    {
                        response  = api.createUser(newusermap.get("username").toString(),newusermap.get("password").toString());
                    }
                    else
                    {
                        //3rd Party OAuth?
                    }
                    
                    break;
                
                    
                /*
                case getAssociatedPerson:
                    
                    response  = api.getAllPerson(data.toString());
                    break;
                */    
                    
                case logout:
                    String key = SecurityUtils.getSubject().getPrincipal().toString();                    
                    response = api.logout();
                    authenticatedMinds.remove(key);
                    break;
                case changePassword:
                    api.changePassword(data.toString());
                    break;
                case learn:
                    api.learn(data);
                    break;
                case log:
                    api.log(data.toString());
                    break;
                case say:
                    api.say(data.toString(), ServerSideAPI.Severity.SUCCESS);
                    break;
                case note:
                    api.note(data.toString());
                    break;
                case alert:
                    api.alert(data.toString());
                    break;
                case convert:
                    response = api.convert(data);
                    break;
                case get:
                    response = api.get(data);
                    if (response == null) response = Void.TYPE;
                    break;
                case set:
                    response = api.set(JSON.mappify(data));
                    if (response == null) response = Void.TYPE;
                    break;
                case create:
                    response = api.create(data);
                    if (response == null) response = Void.TYPE;
                    break;
                case destroy:
                    //XXX: destroy should return status indicator
                    response = api.destroy(data);
                    if (response == null) response = Void.TYPE;
                    break;
                case query:
                    response = api.query(JSON.mappify(data));
                    break;
                case getAll:
                    response = api.getAll(asList(data));
                    break;
                case createAll:
                    response = api.createAll(asList(data));
                    break;
                case destroyAll:
                    api.destroyAll(asList(data));
                    response = Void.TYPE;
                    break;
                case setAll:
                    api.setAll(asList(data));
                    response = Void.TYPE;
                    break;
                case queryOne:
                    response = api.queryOne(JSON.mappify(data));
                    break;

                case run:
                    response = api.run(data);
                    break;
                case listen:
                    api.listen(data.toString());
                    break;
                case unlisten:
                    api.unlisten(data.toString());
                    break;
                case validate:
                    response = api.validate(JSON.mappify(data));
                    break;
                case grant:
                    Map dataMap = JSON.mappify(data);
                    api.grant(
                            JSON.convert(dataMap.get("id"), ObjectId.class),
                            JSON.convert(dataMap.get("user"), String.class),
                            (Set<String>) JSON.convert( dataMap.get("add"), new TypeReference<Set<String>>(){}),
                            (Set<String>) JSON.convert(dataMap.get("remove"), new TypeReference<Set<String>>(){}) 
                            );
                    response = Void.TYPE;
                    break;
                    
                case grantAll:
                    dataMap = JSON.mappify(data);
                    api.grantAll(                            
                            (Set<ObjectId>) JSON.convert(dataMap.get("id"), new TypeReference<Set<ObjectId>>(){}),
                            JSON.convert(dataMap.get("user"), String.class),
                            (Set<String>) JSON.convert( dataMap.get("add"), new TypeReference<Set<String>>(){}),
                            (Set<String>) JSON.convert(dataMap.get("remove"), new TypeReference<Set<String>>(){}) 
                            );
                    response = Void.TYPE;
                    break;
                default:
                    log.warn("Unimplemented channel {}", request.getChannel());
                    response = Void.TYPE;
                    break;
            }
            
            if (response == Void.TYPE)
                connection.deregister(
                    request.getChannel(),
                    request.getRequestId()
                );
            else
                connection.send(new Message(
                    request.getChannel(),
                    response,
                    request.getRequestId(),
                    null
                ));
            
        } catch (Exception e) {
            api.say(e.getMessage(), ServerSideAPI.Severity.FAILURE, request.getRequestId());
            log.error(e.getMessage(), e);
            //deregister failed call
                        //TODO: message client with failure
                connection.deregister(
                    request.getChannel(),
                    request.getRequestId()
                );
        } finally {
            if (ClothoRealm.ANONYMOUS_USER.equals(SecurityUtils.getSubject().getPrincipal())){
                SecurityUtils.getSubject().logout();
            }
        }
    }

    //Start JCA's hack of a pubsub, to be replaced by Ernst
   /* void publish(Sharable object) {
     try {
     System.out.println("Ernst, this needs to be implemented.  Push object via pubsub.");
     String uuid = object.getUUID().toString();
     Map<String, Object> msg = persistor.
     HashSet<WeakReference<ClientConnection>> targets = pubsub.get(uuid);
     for(WeakReference<ClientConnection> wr : targets) {
     ClientConnection conn = wr.get();
     if(conn==null) {
     continue;
     }
     try {
     sendMessage(conn, msg);
     } catch(Exception err) { }
     }
                        
                
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        
     }*/
    void publish(Map object) {
        //doesn't do anything.
        //Is here because Stephanie thinks we should publish the JSON data, not the actual object
    }
    private HashMap<String, HashSet<WeakReference<ClientConnection>>> pubsub = new HashMap<>();

    void register(ClientConnection connection, Sharable object) {
        String uuid = object.getId().toString();
        HashSet<WeakReference<ClientConnection>> existing = pubsub.get(uuid);
        if (existing == null) {
            existing = new HashSet<>();
        }

        existing.add(new WeakReference<>(connection));
        pubsub.put(uuid, existing);
    }
    //End JCA's hack of a pubsub, to be replaced by Ernst
    
    private Mind getMind(ClientConnection connection) {
        String id = connection.getId();
        if (minds.containsKey(id)) {
            Mind mind = minds.get(id);
            if (mind.getConnection() != connection){
                //XXX: this is probably disasterous in some edge cases
                //because jetty preserves the sesson id across websocket close/open, need to check to see if the connection object in the mind is stale
                mind.setConnection(connection);
            }
            return mind;
        }

        Mind mind = new Mind();

        mind.setConnection(connection);
        
        minds.put(id, mind);
        return mind;
    }
    
    private Map<String, Mind> minds;
    private Map<String, Mind> authenticatedMinds = new HashMap<>();

    private Mind getAuthenticatedMind(String username, ClientConnection connection, Message request) {
        //XXX: this whole method is janky
        if (authenticatedMinds.containsKey(username)) {
            return authenticatedMinds.get(username);
        }


        Map<String, Object> query = new HashMap();
        query.put("username", username);
        query.put("schema", Mind.class.getCanonicalName());
        try {
            Iterable<ObjBase> minds = persistor.find(query);

            Mind mind;

            if (!minds.iterator().hasNext()) {
                mind = new Mind();
                authenticatedMinds.put(username, mind);
            } else {
                mind = (Mind) minds.iterator().next();
            }

            return mind;
        } catch (Exception ex) {
            Mind mind = new Mind();
            mind.setConnection(connection);
            authenticatedMinds.put(username, mind);
            
            //tell user mind retrieval failed
            ServerSideAPI api = new ServerSideAPI(mind, persistor, this, realm, request.getRequestId());
            api.say("Mind retrieval encountered an exception: " + ex.getMessage(), ServerSideAPI.Severity.WARNING);
            return mind;
        }
    }
    
    
    //TODO: make everything in the API use iterators instead of lists
    private static List asList(Object o){
        //iterator case
        if (Iterator.class.isInstance(o)){
            return Lists.newArrayList((Iterator) o);
        } else if (Iterable.class.isInstance(o)){
            return Lists.newArrayList((Iterable) o);
        } else if (o instanceof Object[]){
            return Arrays.asList((Object []) o);
        }
        throw new UnsupportedOperationException("Couldn't convert argument to a List");
    }

}
