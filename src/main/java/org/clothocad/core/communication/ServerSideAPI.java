//clotho.show('CT-sample-view', '['static-admin-instance-is-uuid']', '{}'); clotho.show('CT-sample-view', '[]', '{}');
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
 ENHANCEMENTS, OR MODIFICATIONS.
 */
package org.clothocad.core.communication;

import com.fasterxml.jackson.core.JsonParseException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.persistence.EntityExistsException;
import javax.persistence.EntityNotFoundException;
import javax.script.ScriptException;
import javax.validation.ConstraintViolation;
import javax.validation.ConstraintViolationException;
import lombok.Getter;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.authc.AuthenticationException;
import org.apache.shiro.authc.UsernamePasswordToken;
import org.apache.shiro.authz.UnauthorizedException;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.subject.Subject;
import org.clothocad.core.ReservedFieldNames;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.aspects.Interpreter.Interpreter;
import org.clothocad.core.communication.mind.Widget;
import org.clothocad.core.datums.Argument;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Module;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.execution.ConverterFunction;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.execution.ScriptAPI;
import org.clothocad.core.execution.subprocess.SubprocessExec;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.ReflectionUtils;
import org.clothocad.core.security.ClothoAction;
import org.clothocad.core.security.ClothoPermission;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.util.JSON;
import org.clothocad.core.util.XMLParser;
import org.clothocad.model.Person;
import org.python.objectweb.asm.Type;

/**
 * The ServerSideAPI relays the server methods that can be invoked by a client
 * on the server to the Aspects.
 *
 * It gets injected into scripting engines to handle conversion of language. It
 * goes into both the Mind's scripting engine and those of the Executor.
 *
 * It also is for securing the non-secure methods in the Aspects.
 *
 * So, it's just a wrapper really. It's all expressed in terms of Map<String,
 * Object>s and Lists to facilitate conversions. Each different context (Mind,
 * REST, websocket, java client, etc.) is expecting slightly different object
 * representations and synchronization models, so there is necessarily
 * interpretor logic in Router/ Communicator that handles this.
 *
 * API methods should return their result, instead of sending the result as a
 * side effect
 *
 * @author John Christopher Anderson
 */
@Slf4j
public class ServerSideAPI {

    private final Router router;
    private final ClothoRealm realm;
    @Getter
    private final Persistor persistor;
    private final String requestId;
    @Getter
    private final Mind mind;
    private final MessageOptions options;
    private static final Set<Language> executeExternal = EnumSet.of(Language.PYTHON);

    public ServerSideAPI(Mind mind, Persistor persistor, Router router, ClothoRealm realm, String requestId) {
        this(mind, persistor, router, realm, requestId, new MessageOptions());
    }

    public ServerSideAPI(Mind mind, Persistor persistor, Router router, ClothoRealm realm, String requestId, MessageOptions options) {
        this.persistor = persistor;
        this.mind = mind;
        this.requestId = requestId;
        this.router = router;
        this.options = options;
        this.realm = realm;
    }

    public final Iterable<Map<String, Object>> autocomplete(String userText) {
        //This is needed because the subString is in the format {query=[subString]}
        userText = userText.substring(7, userText.length() - 1);

        //Add the word suggestions from the global Trie
        return persistor.getCompletions(userText);
    }
    //JCA:  works pushing a dummy message to the client, probably should be wrapped into get(...)

    public final String autocompleteDetail(String uuid) {
        try {
            Map<String, Object> msg = JSON.deserializeObjectToMap("{\"channel\":\"autocompleteDetail\",\"data\":{\"uuid\":\"1234567890\",\"text\":\"This is a command\",\"command\":\"clotho.run('230sdv-232', '18919e-18')\",\"versions\":[{\"uuid\":\"uuid123\",\"text\":\"Reverse Complement Tool\",\"author\":{\"uuid\":\"uuid_author_123\",\"name\":\"Joe Schmo\",\"email\":\"joe@schmo.com\",\"biography\":\"This is a biography about Joe Schmo. It's not too long. \"},\"description\":\"Aenean lacinia bibendum nulla sed consectetur. Cum sociis natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Donec ullamcorper nulla non metus auctor fringilla. Maecenas faucibus mollis interdum. Etiam porta sem malesuada magna mollis euismod.\",\"usage\":{\"executed\":\"35\",\"successful\":\"27\",\"positive\":\"12\",\"negative\":\"3\"}},{\"uuid\":\"uuid456\",\"text\":\"pBca 1256\",\"author\":{\"uuid\":\"uuid_author_456\",\"name\":\"Chris Anderson\",\"email\":\"chris@anderson.com\",\"biography\":\"This is a biography about Chris Anderson. It's different than Joe's... It's a little longer. Yada yada yada. Here's some latin. It should get truncated on the server or we could write our own directive to handle truncating (easy). Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.\"},\"description\":\"Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.\",\"usage\":{\"executed\":\"8\",\"successful\":\"8\",\"positive\":\"6\",\"negative\":\"0\"}}]}}");
            return msg.get("data").toString();
        } catch (JsonParseException ex) {
            ex.printStackTrace();
            return null;
        }

    }

    //clotho.run("aa7f191e810c19729de86101", ["53581f9e9e7d7a2fda8c36a7"]);   revcomp pBca1256
    //clotho.run("aa7f191e810c19729de86101", ["atcg"]);  revcomp atcg
    public final Object submit(Object data) {
        //Extract the query String and tokens
        Map<String, Object> json = JSON.mappify(data);
        String query = (String) json.get("query");
//        List<Object> clientTokens = (List<Object>) json.get("tokens");

//        //Resolve the commands to tokens
//        System.out.println("++ The command submitted is: " + query);
//        String[] tokens = query.split("\\s+");
//        for(String str : tokens) {
//            System.out.println("token: " + str);
//        }
//
//        Object out = null;
//        try {
//            return tryRun(tokens);
//        } catch(Exception err) {
//        }
//        
//        try {
//            return trySingleWord(tokens);
//        } catch(Exception err) {
//        }
//        
//        try {
//            return tryAPIWord(tokens);
//        } catch(Exception err) {
//        }
//        
        //Run the command assuming it's javascript
        //say(command, Severity.MUTED, null, true);
        try {
            Object returnValue = mind.runCommand(query, getScriptAPI());
            //If the command successfully executed, it gets retained
            mind.addLastCommand(Channel.submit, query);
            return returnValue;
        } catch (ScriptException ex) {
            //disambiguate(command);  //JCA:  temporarily disabled for testing, also not fully hooked up
            logAndSayError("Unable to process the request", ex);
            return Void.TYPE;
        }
    }

    /**
     * Interprets tokens of a submit as a clotho.run call where the first token
     * is the Function, and the others are args
     *
     * @param tokens
     * @return the result or null if it failed to execute
     */
//    private Object tryRun(String[] tokens) throws Exception {
//        System.out.println("+++  try RUN on args");
//        Function function = null;
//        List<Object> args = new ArrayList<Object>();
//        
//        //Assume the first token is a Function call
//        List<Map> completions = persistor.getCompletions(tokens[0]);
//        String uuid = (String) completions.get(0).get("id").toString();
//        function = persistor.get(Function.class, new ObjectId(uuid));
//        if(function==null) {
//            throw new Exception();
//        }
//        
//        //Iterate through each subsequent token and convert to object if required
//        for(int i=1; i<tokens.length; i++) {
//            try {
//                Object obj = resolveSloppy(tokens[i]);
//                args.add(obj);
//            } catch(Exception err) {
//                args.add(tokens[i]);
//            }
//        }
//        
//        
////JCA TODO: remove this section when converters become available
//        //Extract sequence strings when an object is given as argument for a String field
//        Argument[] funcargs = function.getArgs();
//        for(int i=0; i<funcargs.length; i++) {
//            try {
//                Argument arg = funcargs[i];
//                String str = arg.getType().toString();
//                if(str.equals("class java.lang.String")) {
//                    Map map = (Map) args.get(i);
//                    String newvalue = (String) map.get("sequence");
//                    args.remove(i);
//                    args.add(i, newvalue);
//                }
//            } catch(Exception err) {
//                
//            }
//        }
//
//        //Invoke run with the ScriptEngine
//        System.out.println("Running the resolved sloppy arguments with function:\n " + function.toString());
//        System.out.println("And args: " );
//        for(Object obj : args) {
//            System.out.println(obj.toString());
//        }
//        Object out = null;
//        out = run(function, args);
//        return out;
//    }
//    
//    private Object trySingleWord(String[] tokens) {
//        System.out.println("+++  try Single Word get");
//        return resolveSloppy(tokens[0]);
//    }
//    
//    private Object resolveSloppy(String word) {
//        List<Map> completions = persistor.getCompletions(word);
//        String uuid = (String) completions.get(0).get("id");
//        return get(uuid);
//    }
//    
//    private Object tryAPIWord(String[] tokens) throws Exception {
//        System.out.println("+++  try first word is API word");
//        String firstWord = tokens[0].toLowerCase();
//        
//        //Example:  get pBca1256
//        if(firstWord.equals("get")) {
//            //Interpret the next token as a clotho.get request
//            return resolveSloppy(tokens[1]);
//            
//        } else if(firstWord.equals("run")) {
//            String[] newtok = new String[tokens.length-1];
//            for(int i=0; i<newtok.length; i++) {
//                newtok[i] = tokens[i+1];
//            }
//            return tryRun(newtok);
//        } else {
//            throw new Exception();
//        }
//    }
//    
    public final void learn(Object data) {
        //might already be data?
        Map<String, Object> json = JSON.mappify(data);
        learn((String) json.get("userInput"), (Message) json.get("command"));
    }

    //clotho.learn("test1", "clotho.say('Hi there');")
    public final void learn(String userInput, Message command) {

        if (command == null) {
            command = mind.getLastCommands().get(0);
            if (command == null) {
                System.out.println("This is something that will happen if you haven't issued a previous command, it has nothing to learn");
                return;
            }
            Interpreter.get().learnNative(userInput, JSON.mappify(command));
        }
    }

    /*needs to be replaced by general-purpose createAccount method, see account
     creation discussion in email thread on dev list 
      
     currently assumes that new person should be created with no data - generic 
     method should optionally accept new person data*/
    public final Map<String, Object> createUser(String username, String password) {

        Subject subject = SecurityUtils.getSubject();
        if (!ClothoRealm.ANONYMOUS_USER.equals(subject.getPrincipal())) {
            say("You are already logged in as " + subject.getPrincipal() + "!", Severity.FAILURE);
            return null;
        }
        try {
            realm.addAccount(username, password);
            say("New user " + username + " created.", Severity.SUCCESS);

            subject.logout();
            subject.login(new UsernamePasswordToken(username, password));

            Person newPerson = new Person(username);
            newPerson.setPrimaryEmail(username);
            newPerson.setEmailAddress(username);
            newPerson.setId(new ObjectId());
            ObjectId personId = persistor.save(newPerson);

            realm.addPrincipal(username, personId, "persons");

            Map<String, Object> result = new HashMap<>();
            result.put("id", personId);
            result.put("accessToken", "dummy");
            result.put("app_id", "dummy");
            return result;
        } catch (EntityExistsException e) {
            say("User " + username + " exists.", Severity.FAILURE);
            return null;
        } catch (Exception e) {
            logAndSayError("There was a problem creating the new account.", e);
            return null;
        } finally {
            subject.logout();
            subject.login(ClothoRealm.getAnonymousUserToken());
        }
    }

    public final Map<String, Object> login(String username, String password) {
        ObjectId userId = null;
        Subject subject = SecurityUtils.getSubject();

        if (!subject.isAuthenticated()) {
            try {
                subject.login(new UsernamePasswordToken(username, password));
            } catch (AuthenticationException e) {
                logAndSayError("Authentication attempt failed for username " + username, e);
                return null;
            }

            say("Welcome, " + username, Severity.SUCCESS);
            Map<String, Object> result = new HashMap<>();
            Object personId;
            Collection personPrincipals = subject.getPrincipals().fromRealm("persons");
            if (personPrincipals.isEmpty()) {
                personId = null;
            } else {
                personId = personPrincipals.iterator().next();
            }
            result.put("id", personId);
            result.put("accessToken", "dummy");
            result.put("app_id", "dummy");
            log.info("User {} logged in", username);

            return result;
        } else {
            say("Please logout first.", Severity.FAILURE);
            return null;
        }

    }

    public final Object loginOAuth(String username, Object Map) {
        ObjectId userId = null;
        return userId;
    }

    public final boolean logout() {
        if (SecurityUtils.getSubject().isAuthenticated()) {
            String username = SecurityUtils.getSubject().getPrincipal().toString();
            mind.setUsername(username);
            //XXX: need some kind of error recovery if mind save fails
            try {
                persistor.save(mind);
            } catch (Exception e) {
                say("There was a problem saving your mind. You will still be logged out, but some settings may not be saved.", Severity.WARNING);
            }
            SecurityUtils.getSubject().logout();
            say("Logged out", Severity.SUCCESS);
            log.info("User {} logged out", username);
            return true;
        }

        say("You are not logged in", Severity.WARNING);
        return false;
    }

    public final boolean changePassword(String newPassword) {
        return true;
    }

    public final void clear() {
        mind.clear();
        say("The mind has been cleared", Severity.SUCCESS);
    }

    /**
     *
     * @param message
     * @param severity "text-error", "text", "text-warning", "text-success" see
     * search-directives.js for 'from', is client or server
     */
    public final void say(String message, Severity severity) {
        say(message, severity, null, false);

    }

    public final void say(String message, Severity severity, String recipients) {
        say(message, severity, recipients, false);

    }

    protected void say(String message, Severity severity, String recipients, boolean isUser) {
        //if say is turned off in options, send nothing
        if (options.isMute()) {
            return;
        }

        //Resolve the recipients
        //XXX: doesn't currently handle multiple recipients
        //List<Sharable> listUsers = resolveToExistentSharablesList(recipients);
        String source;
        if (isUser) {
            source = "client";
        } else {
            source = "server";
        }
        Map<String, Object> data = new HashMap<>();
        data.put("text", message);
        data.put("from", source);
        data.put("class", severity);
        data.put("timestamp", new Date().getTime());

        Message msg = new Message(Channel.say, data, requestId, null);
        router.sendMessage(mind.getConnection(), msg);
    }

    //JCA:  Java side looks fine, but client code crashes browser
    //clotho.alert("this is an alert!");
    public final void alert(String message) {
        router.sendMessage(
                mind.getConnection(),
                new Message(Channel.alert, message, null, null));
    }

    //JCA:  This runs, and the message goes to the console.log spot.
    //clotho.log("I did some minipreps today");
    public final void log(String message) {
        log.debug("log has: {}", message);
        router.sendMessage(
                mind.getConnection(),
                new Message(Channel.log, message, null, null));
    }

    //Make note of this message in my notebook
    public final void note(String message) {
        System.out.println("I need to put this in your notebook, but i'm not implemented");
        //These have the same structure as a say, but they are stored.

        say("I've stored your note (but not really): " + message, Severity.MUTED);
    }

    protected void send(Message message) {
        router.sendMessage(mind.getConnection(), message);
    }

    private Object unwrap(Object o) {
        while (o instanceof Collection) {
            Iterator iterator = ((Collection) o).iterator();
            if (iterator.hasNext()) {
                o = iterator.next();
            } else {
                throw new IllegalArgumentException("Cannot unwrap an empty collection");
            }
        }

        return o;
    }

    public Map<String, Object> get(Object o) {
        o = unwrap(o);
        return get(new ObjectId(o));
    }

    public Map<String, Object> get(ObjectId id) {
        try {
            Map<String, Object> out = persistor.getAsJSON(id, options.getPropertiesFilter());
            say(String.format("Retrieved object #%s", id.toString()), Severity.SUCCESS);
            return out;

        } catch (UnauthorizedException e) {
            say(String.format("The current user does not have read access for object #%s", id.toString()), Severity.FAILURE);
            return null;
        } catch (EntityNotFoundException e) {
            say(String.format("No object with id %s", id.toString()), Severity.FAILURE);
            return null;
        }
    }

    public final List<Map<String, Object>> getAll(List objects) {
        List<Map<String, Object>> returnData = new ArrayList<>();
        for (Object obj : objects) {
            try{
                returnData.add(get(obj));
            } catch (Exception e){
                returnData.add(null);
            }
        }
        return returnData;
    }

    private ObjectId resolveId(String id) {
        ObjectId uuid;
        try {
            uuid = new ObjectId(id);
        } catch (IllegalArgumentException e) {
            say(e.getMessage(), Severity.FAILURE);
            return null;
        }
        return uuid;
    }

    private void publish(Map obj) {
        router.publish(obj);
    }

    public final List<ObjectId> setAll(List<Map<String, Object>> values) {
        List<ObjectId> out = new ArrayList<>();
        for (Map<String, Object> spec : values) {
            out.add(set(spec));
        }
        return out;
    }

    public final ObjectId set(Map<String, Object> values) {
        try {

            if (values.get("id") == null) {
                say("set: No uuid provided", Severity.WARNING);
                return create(values);
            }

            ObjectId uuid = resolveId(values.get("id").toString());
            if (uuid == null) {
                return null;
            }

            if (!persistor.has(uuid)) {
                say("set: No object with this id exists", Severity.WARNING);
                return create(values);
            }

            //Grab the object to be altered
            Map<String, Object> original = persistor.getAsJSON(uuid);

            //TODO: run validation check on saved values
            persistor.save(values);

            //Confirm that the new data is different than the old data
            Map<String, Object> altered = persistor.getAsJSON(uuid);
            if (original.toString().equals(altered.toString())) {
                say(String.format("Object #%s named %s was unmodified.", uuid.toString(), altered.get("name")), Severity.WARNING);
            }

            //Contact the user to notify them that they modified an object
            say(String.format("Successfully modified object #%s named %s", uuid.toString(), altered.get("name")), Severity.SUCCESS);

            //Relay the data change to listening clients
            publish(altered); //publish by uuid?
            return uuid;
        } catch (UnauthorizedException e) {
            say(e.getMessage(), Severity.FAILURE);
            return null;
        } catch (Exception e) {
            logAndSayError(String.format("Error setting %s: %s", values.toString(), e.getMessage()), e);
            return null;
        }
    }

    public final String create(Object o) {
        ObjectId id = create(JSON.mappify(o));
        if (id == null) {
            return null;
        }
        return id.toString();
    }

    //TODO: some global solution for jsonifying ObjectIds
    public final List<String> createAll(List<Object> objects) {
        List<String> returnData = new ArrayList<>();
        //list of selectors?
        for (Object obj : objects) {
            try{
                ObjectId id = create(JSON.mappify(obj));
                returnData.add(id == null? null : id.toString());               
            } catch (Exception e){
                returnData.add(null);
            }
        }
        return returnData;
    }

    public ObjectId create(Map<String, Object> obj) {

        try {
            //Confirm that there is no pre-existing object with this uuid
            String idKey = null;
            if (obj.containsKey("id")) {
                idKey = "id";
            }

            if (idKey != null) {

                ObjectId uuid = resolveId(obj.get(idKey).toString());
                if (uuid == null) {
                    return null;
                }
                if (persistor.has(uuid)) {
                    say("An object with the uuid " + uuid + " already exists.  No object was created.", Severity.FAILURE);
                    return null;
                }

                obj.put(idKey, new ObjectId(obj.get(idKey).toString()));
            }


            try {
                ObjectId id = persistor.save(obj);
                //TODO: Relay the data change to listening clients
                //TODO: create sets author to current user

                //Return the JSON of the new object as a String
                say(String.format("Created object #%s named %s", id.toString(), obj.get("name")), Severity.SUCCESS);
                return id;
            } catch (ConstraintViolationException e) {
                say(String.format("Validation failed: %s. No object was created.", e.getMessage()), Severity.FAILURE);
                return null;
            }

        } catch (UnauthorizedException e) {
            say("The current user does not have write access for this domain", Severity.FAILURE);
            return null;
        }

        //catch (Exception e) {
        //    logAndSayError(String.format("Error creating %s: %s", obj.toString(), e.getMessage()), e);
        //    return null;
        //}
    }

    protected void logAndSayError(String message, Exception e) {
        log.error(message, e);
        say(message, Severity.FAILURE);
    }

    public final ObjectId destroy(Object id) {
        if (id == null) {
            return null;
        }

        ObjectId resolvedId = new ObjectId(id);
        try {
            try {
                persistor.delete(resolvedId);
            } catch (UnauthorizedException e) {
                say(e.getMessage(), Severity.FAILURE);
                return null;
            }
            say(String.format("Destroyed object #%s", resolvedId.toString()), Severity.SUCCESS);
            return resolvedId;
        } catch (Exception e) {
            logAndSayError(String.format("Error destroying %s: %s", id.toString(), e.getMessage()), e);
            return null;
        }
    }

    public final List<ObjectId> destroyAll(List<Object> objects) {
        List<ObjectId> out = new ArrayList<>();
        for (Object obj : objects) {
            out.add(destroy(obj));
        }
        return out;
    }

    public List<Map<String, Object>> query(Map<String, Object> spec) {
        List<Map<String, Object>> objs;
        try {
            //Relay the query to Persistor and return the hits
            objs = persistor.findAsJSON(spec, options.getPropertiesFilter(), options.getMaxResults());
            say("Found " + objs.size() + " objects that satisfy your query", Severity.SUCCESS);
            return objs;
        } catch (Exception e) {
            logAndSayError(String.format("Error querying %s: %s", spec.toString(), e.getMessage()), e);
            e.printStackTrace();
            return new ArrayList<>();
        }
    }

    public final Object convert(Object o) {
        final Map<String, Object> data = JSON.mappify(o);
        Object result = null;
        Object convertThis = null;
        Schema targetSchema = null;
        Schema currentSchema = null;
        for (Map.Entry<String, Object> entry : data.entrySet()) {
            if (entry.getKey().equals("convertTo")) {
                targetSchema = persistor.get(Schema.class, new ObjectId(((HashMap) (entry.getValue())).get("id")));
            } else if (entry.getKey().equals("convert")) {
                convertThis = entry.getValue();
            }
        }
        if (convertThis == null || targetSchema == null) {
            if (convertThis == null) {
                say("Object to convert to was null", Severity.FAILURE);
            }
            if (targetSchema == null) {
                say("Convert To Schema was null", Severity.FAILURE);
            }
        }
        {
            result = convert(convertThis, targetSchema);
        }
        return result;
    }

    public final Object convert(Object obj, Schema schema) {
        Object result = null;
        String idVal = (String) ((HashMap) obj).get("id");
        Map<String, Object> objMap = persistor.getAsJSON(new ObjectId(idVal));
        Schema currentSchema = persistor.get(Schema.class, new ObjectId(objMap.get("schema")));


        boolean foundFunc = false;
        ConverterFunction resultCFunc = null;
        Collection<ConverterFunction> convlist = persistor.getAll(ConverterFunction.class);
        for (ConverterFunction xconvfunc : convlist) {
            if (xconvfunc.convertTo.equals(schema)) {
                if (xconvfunc.convertFrom.equals(currentSchema)) {
                    foundFunc = true;
                    resultCFunc = xconvfunc;
                    String convfuncName = xconvfunc.getName();
                    if (convfuncName != null) {
                        say("Converter Funcion : " + xconvfunc.getName() + "found and will now be executed.", Severity.SUCCESS);
                    } else {
                        say("Converter Function Found, but no name found", Severity.WARNING);
                    }
                    break;
                }
            }
        }
        List<Object> args = new ArrayList<Object>();
        args.add(obj);
        if (foundFunc) {
            if (resultCFunc.getFunction() == null) {
                say("Converter has a null function.", Severity.FAILURE);
            } else {
                try {
                    result = run(resultCFunc.getFunction(), args);
                } catch (ScriptException ex) {
                    log.error("", ex);
                }
            }
        } else {
            say("No suitable Converter Found.", Severity.FAILURE);
        }


        return result;
    }

    private Object runAsSubprocess(final Function function, final List<Object> args) {
        final Map<String, Object> funcJSON =
                persistor.getAsJSON(function.getId());
        Object out = SubprocessExec.run(
                this,
                funcJSON,
                args,
                new SubprocessExec.EventHandler() {
            @Override
            public void onFail(final byte[] err) {
                say_helper(err, Severity.FAILURE);
            }

            @Override
            public void onSuccess(final byte[] err) {
                say_helper(err, Severity.NORMAL);
            }

            private void say_helper(final byte[] err,
                    final Severity sev) {
                if (err.length != 0) {
                    say(new String(err, StandardCharsets.UTF_8), sev);
                }
            }
        });
        return out;
    }

    public List<Object> resolveIds(Argument[] functionArguments, List<Object> suppliedArguments) {
        //Function spec might be incomplete - might have not declared argument types at all;
        if (functionArguments == null) {
            return suppliedArguments;
        }

        suppliedArguments = new ArrayList(suppliedArguments); //because some list implementations don't support set
        for (int i = 0; i < functionArguments.length && i < suppliedArguments.size(); i++) {
            //XXX: can't actually declare polymorphic types yet - simple case only
            Class declaredType = functionArguments[i].getType();
            Object argument = suppliedArguments.get(i);
            if (declaredType != null && ObjBase.class.isAssignableFrom(declaredType) && argument != null) {
                ObjectId idArgument;

                if (String.class.isInstance(argument)) {
                    idArgument = new ObjectId(argument);
                } else if (ObjectId.class.isInstance(argument)) {
                    idArgument = (ObjectId) argument;
                } else {
                    continue;
                }
                if (persistor.has(idArgument)) {
                    suppliedArguments.set(i, persistor.getAsJSON(new ObjectId(argument)));
                }
            }
        }

        return suppliedArguments;
    }

    //TODO: needs serious cleaning up
    public final Object run(Object o)
            throws ScriptException,
            IllegalAccessException,
            IllegalArgumentException,
            InvocationTargetException {

        final Map<String, Object> data = JSON.mappify(o);
        final List<Object> args;

        if (data.get(ID) == null) {
            say("No runnable specified (id missing)", Severity.FAILURE);
            return Type.VOID;
        }

        //ensure args is list & exists
        try {
            args = (List) data.get("args");
        } catch (ClassCastException e) {
            say("Arguments must be a list or array", Severity.FAILURE);
            return Type.VOID;
        }

        //XXX: check that functionData exists
        Map<String, Object> functionData = persistor.getAsJSON(new ObjectId(data.get(ID)), null, true);


        if (functionData.containsKey("schema") && (functionData.get("schema").toString().endsWith("Function")
                || functionData.get("schema").toString().endsWith("Module"))) {

            Module module = persistor.get(Module.class, new ObjectId(data.get(ID)), true);

            if (executeExternal.contains(module.getLanguage())) {
                //execute using process launcher
                try {
                    Function function = (Function) module;
                    return runAsSubprocess(function, resolveIds(function.getArgs(), args));
                } catch (ClassCastException e) {
                    logAndSayError("Can only execute Functions in subprocess: ", e);
                    return Void.TYPE;
                }

            } else {
                //execute using script engine
                try {
                    if (data.get("function") != null) {
                        // was a function to execute indicated? If so, execute as module w/ the function named as target
                        String functionName = data.get("function").toString();
                        Function function = module.getFunction(functionName);
                        return mind.invokeMethod(module, functionName, resolveIds(function.getArgs(), args), getScriptAPI());

                    } else { // if not, execute as Function
                        try {
                            Function function = (Function) module;
                            return mind.invoke(function, resolveIds(function.getArgs(), args), getScriptAPI());
                        } catch (ClassCastException e) {
                            logAndSayError(data.get(ID).toString() + " is not itself a Function, but no method name provided in 'function' field", e);
                            return Void.TYPE;
                        }
                    }

                } catch (ScriptException e) {
                    logAndSayError("Script Exception thrown: " + e.getMessage(), e);
                    return Void.TYPE;
                } catch (NoSuchMethodException ex) {
                    logAndSayError("No such function found", ex);
                    return Void.TYPE;
                }
            }
        }

        //"Old" run, for non-Module objects

        //resolve any ids to their objects
        if (args != null) for (int i = 0; i < args.size(); i++) {
            try {
                ObjectId id = new ObjectId(args.get(i).toString());
                //must have read privs on args
                args.set(i, persistor.get(ObjBase.class, id));
            } catch (EntityNotFoundException e) {
                //XXX: warn here? fail here?
            }
        }

        //check for permissions
        ObjectId id = new ObjectId(data.get(ID).toString());
        persistor.checkPriv(id, "run");
        //reflectively (ugh) run function of instance
        ObjBase instance = persistor.get(ObjBase.class, id, true);


        Method method = ReflectionUtils.findMethodNamed(data.get("function").toString(), args.size(), instance.getClass());
        Object result = method.invoke(instance, args.toArray());
        if (method.getReturnType().equals(Void.TYPE)) {
            return Void.TYPE;
        }
        List<Object> results = new ArrayList<>();
        if (result instanceof Iterable) {
            for (Object r : ((Iterable) result)) {
                if (r instanceof ObjBase) {
                    results.add(persistor.save((ObjBase) r));
                } else {
                    results.add(r);
                }
            }
        } else {
            if (result instanceof ObjBase) {
                return persistor.save((ObjBase) result);
            } else {
                return result;
            }
        }
        return results;
    }

    public final Object run(Function function, List<Object> args) throws ScriptException {
        if (function.getLanguage().equals(Language.PYTHON)) {
            return runAsSubprocess(function, args);
        }

        return mind.evalFunction(function.getCode(), function.getName(), args, getScriptAPI());
    }

    /**
     * Relay method for receiving a show command
     *
     * @param viewRef
     * @param sharables
     * @param sloppy_args
     */
    public final void show(String viewIod, String sharables, String position, String recipients) {
        System.out.println("Show is called " + position);
        /**
         * *
         * try {
         *
         * //IF THE PROCESS WAS TIED TO A COMMAND-INITIATED PROCESS THERE WOULD
         * BE A PARENT DOO THAT NEEDS TO BE FOUND. Doo parentDoo =
         * Hopper.get().extract(null); ShowDoo doo = new ShowDoo(parentDoo);
         *
         * //Gather up referenced objects View view = (View)
         * resolveToObjBase(viewId); List<Sharable> shareList =
         * resolveToSharables(sharables);
         *
         * //Create the widget and put it on its page Page targetPage =
         * mind.getConfig().getPage(socket_id); Widget widget = new
         * Widget(targetPage, view); targetPage.addWidget(widget);
         *
         * //Just for record keeping doo.viewId = view.getId(); doo.widgetId =
         * widget.getId(); doo.collectSharables(shareList);
         *
         * //Create all the commands for the client doo.commandMessageArray =
         * new ArrayList(); doo.commandMessageArray.put(makeCollect(shareList));
         * doo.commandMessageArray.put(makeShowWidget(widget, position));
         * doo.commandMessageArray.put(makeUpdate(widget, shareList,
         * socket_id)); doo.commandMessageArray.put(makeCallback(doo.getId(),
         * socket_id));
         *
         * //Put the doo into the hopper to await a callback
         * Hopper.get().add(doo);
         *
         * //Send the commands Communicator.get().sendClientMessage(socket_id,
         * SendChannels.commandList, doo.commandMessageArray.toString());
         *
         * //Save everything whose state was changed //JCA: THIS NEEDS A
         * CALLBACK/FAILURE RESPONSE THAT REVERTS THIS (EVENTUALLY)
         * Persistor.get().persistObjBase(mind); } catch (Exception e) {
         * Logger.log(Logger.Level.WARN, "", e); e.printStackTrace(); //REVERT
         * THE MIND //SEND CLIENT MESSAGE TELLING THAT IT FAILED } *
         */
    }

    /**
     * Invoke the hard-coded editor appropriate for this particular object *For
     * an Instance, invoke the view set as the default view for this type of
     * data. That's one option. The other is to build up an editor
     * programmatically from the components. Not sure. Actually, I like this
     * second option best. Maybe. I don't know. Second option is more secure.
     *
     * Both options need to be available, I think.
     *
     * @param sharableRef
     */
    public final void edit(String sharableRef) {
        /**
         * try { Sharable sharable = (Sharable) resolveToObjBase(sharableRef);
         *
         * switch (sharable.type()) { case SCHEMA: Schema schema = (Schema)
         * sharable; //Pop up the page, push in data return; case INSTANCE:
         * return; default: return; } //Pop up a new window of } catch
         * (Exception err) { err.printStackTrace(); }
         *
         */
    }

    public final void listen(String args) {
        say("not yet implemented", Severity.FAILURE);
    }

    public final void unlisten(String data) {
        say("not yet implemented", Severity.FAILURE);
    }

//    public final void test() {
//        try {
//            System.out.println("Test has been invoked");
//            Map<String, Object> trail = new HashMap<String, Object>("{\"uuid\":\"trail_123\",\"title\":\"Biosafety Module\",\"author\":\"UC Berkeley\",\"description\":\"<blockquote><p>This is a module on Biosafety. You'll learn about Cras sit amet nibh libero, in gravida nulla. Nulla vel metus scelerisque ante sollicitudin commodo. Cras purus odio, vestibulum in vulputate at, tempus viverra turpis.</p><p><small>Maxwell Bates</small></p></blockquote>\",\"contents\":[{\"module_title\":\"The Basics\",\"pavers\":[{\"paver_title\":\"Introduction\",\"type\":\"template\",\"template\":\"/app/partials/trail_123_1.html\"}]},{\"module_title\":\"Biosafety Levels\",\"pavers\":[{\"paver_title\":\"Introduction\"},{\"paver_title\":\"Biosafety Level 1-2\"},{\"paver_title\":\"Biosafety Level 3-4\"}]},{\"module_title\":\"Assessment\",\"pavers\":[{\"paver_title\":\"Review\"}],\"assessment\":[{\"type\":\"quiz\",\"title\":\"Final Quiz\"}]}]}");
//            Trail i = new Trail("Biosafety Module", trail);
//
//            System.out.println("Trail i has been created: " + i.toString());
//            
//            Persistor.get().save(i);
//            ObjectId id = i.getUUID();
//            
//            String uuid = id.toString();
//
//            say(id.toString());
//            
//            ObjBase reclaimed =  Collector.get().getObjBase(uuid);
//            
//            System.out.println("After re-retrieval from db I have: " + reclaimed.toString());
//        } catch (JSONException ex) {
//            java.util.logging.Logger.getLogger(ServerSideAPI.class.getName()).log(Level.SEVERE, null, ex);
//        }
//    }
    /**
     * Return the Person object associated with the Mind into which this SS API
     * has been injected
     *
     * @return
     */
    private Person getPerson() {
        return mind.getPerson();
    }

    public final void newPage() throws Exception {
        /**
         * *
         * try { List commandMessageArray = new ArrayList();
         * commandMessageArray.put(makeNewPage(UUID.randomUUID().toString()));
         *
         * // // // //Put the doo into the hopper to await a callback //
         * Hopper.get().add(doo);
         *
         * //Send the commands Communicator.get().sendClientMessage(socket_id,
         * SendChannels.commandList, commandMessageArray.toString());
         *
         * //Save everything whose state was changed //JCA: THIS NEEDS A
         * CALLBACK/FAILURE RESPONSE THAT REVERTS THIS (EVENTUALLY)
         * Persistor.get().persistObjBase(mind); } catch (Exception e) {
         * Logger.log(Logger.Level.WARN, "", e); e.printStackTrace(); //REVERT
         * THE MIND //SEND CLIENT MESSAGE TELLING THAT IT FAILED } *
         */
    }

    /**
     * //if a command has a callback, handle it here //minimally removes the
     * doo from the dooHopper //jca: i don't know that i like this being in the
     * api like this. it really //should be on a more secure channel because
     * this all deals with Doos, which //need to be accurately tracked. if you
     * can directly call this callback from //the api, you can just type it into
     * the commmand bar and lie to clotho, which //is no good. so, yeah, it
     * shouldn't be implemented like this, but we can fix that //in the future
     * when we pick through the client-->server communications standard. //also,
     * this is a mind-agnostic thing -- this is what gets called just to confirm
     * //that a process has completed, and everything should probably call this.
     * The doo //hopper maybe should even be an Aspect with this being relayed
     * there. Yeah. That's //it.
     */
    public final void notify(String dooID) {
        /**
         * System.out.println("calling back!"); //Grab the Doo and perhaps
         * perform any commands embedded in there
         *
         * //Terminate it from the Hopper Hopper.get().terminate(dooID);
         *
         */
    }

    /**
     * Instantiate a new workspace page. Only the Proctor can construct
     * arguments for TRAILS or specific editors. Security thing.
     *
     * @return
     * @throws Exception
     */
    public static Map<String, Object> makeNewPage(String page_id) throws Exception {
        Map<String, Object> out = new HashMap<>();
        out.put("mode", "WORKSPACE");
        out.put("ephemeral_link_page_id", page_id);
        out.put("command", "addPage");
        return out;
    }

    /**
     * Instruct the client to display the GUI of a view and position it
     * according to positioning parameters
     *
     * @param viewId
     * @param sloppy_args
     * @return
     * @throws Exception
     */
    public Map<String, Object> makeShowWidget(Widget widget, String position) throws Exception {
//        Map<String, Object> out = new HashMap<String, Object>();
//        Map<String, Object> positioninfo = resolveToMap<String, Object>(position);
//        View view = widget.getView();
//
//        //Load the html and js scripts into a Map<String, Object>, insert the uuid for the widget
//        String widgetId = widget.getId();
//        System.out.println("graphics script: "+view.getGraphicsScript());
//        String html = replaceWidgetId(view.getGraphicsScript(), widgetId);
//        String onshow = replaceWidgetId(view.getOnShowScript(), widgetId);
//
//        out.put("widget_id", widgetId);
//        out.put("content", html);
//        out.put("on_show", onshow);
//        out.put("parent_widget_id", "widget_space_root");
//        out.put("command", "showWidget");
//
//        return out;
        return null;
    }

    /**
     * Replace the _widget_id phrases with the actual uuid
     *
     * @param script
     * @param widgetIdPrefix
     * @return
     */
    public static final String replaceWidgetId(String script, String widgetIdPrefix) {
        try {
            return XMLParser.addPrefixToTagAttribute(
                    script,
                    ReservedFieldNames.ID,
                    widgetIdPrefix);
        } catch (Exception ex) {
            log.error("", ex);
        }
        return null;
    }

    Object queryOne(Map<String, Object> query) {
        List result = persistor.findAsJSON(query, options.getPropertiesFilter(), 1);
        if (result.isEmpty()) {
            return Void.TYPE;
        }
        return result.get(0);
    }

    Set<ConstraintViolation<?>> validate(Map<String, Object> data) {
        try {
            persistor.validateBSON(data);
        } catch (IllegalArgumentException iae) {
            say(String.format("Could not validate: %s", iae.getMessage()), Severity.WARNING);

        } catch (ConstraintViolationException e) {
            say(String.format("Validation unsuccessful: %s", e.getMessage()), Severity.FAILURE);
            return e.getConstraintViolations();
        }
        say("Validation successful.", Severity.SUCCESS);
        return new HashSet<>();
    }

    public void grantAll(Collection<ObjectId> ids, String principal, Set<String> addPermissions, Set<String> removePermissions) {
        for (ObjectId id : ids) {
            grant(id, principal, addPermissions, removePermissions);
        }
    }

    public void grant(ObjectId id, String principal, Set<String> addPermissions, Set<String> removePermissions) {
        try {
            Set<ClothoAction> add = new HashSet<>();
            Set<ClothoAction> remove = new HashSet<>();
            for (String permString : addPermissions) {
                try {
                    ClothoPermission perm = ClothoPermission.valueOf(permString.toUpperCase());
                    if (perm.equals(ClothoPermission.PUBLIC)) {
                        //add public status
                        realm.setPublic(id);
                    } else {
                        add.addAll(perm.actions);
                    }
                } catch (IllegalArgumentException e) {
                    //wrong permstring name
                    logAndSayError("Invalid permission name", e);
                }
            }
            realm.addPermissions(principal, add, id);
            for (String permString : removePermissions) {
                try {
                    ClothoPermission perm = ClothoPermission.valueOf(permString.toUpperCase());
                    if (perm.equals(ClothoPermission.PUBLIC)) {
                        //add public status
                        realm.removePublic(id);
                    } else {
                        remove.addAll(perm.removedActions);
                    }
                } catch (IllegalArgumentException e) {
                    logAndSayError("Invalid permission name", e);
                }
            }
            realm.removePermissions(principal, remove, id);
        } catch (AuthorizationException e) {
            logAndSayError("Cannot grant/remove permission", e);
        }
    }

    public static enum Severity {

        SUCCESS,
        WARNING,
        FAILURE,
        NORMAL,
        MUTED
    }

    private ScriptAPI getScriptAPI() {
        return new ScriptAPI(this);
    }
}
