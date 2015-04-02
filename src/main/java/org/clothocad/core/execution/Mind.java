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

package org.clothocad.core.execution;

import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.script.ScriptException;
import lombok.Getter;
import lombok.Setter;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Module;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.model.Person;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A 'Mind' is the thing on the server that connects up a User and a Client.  Minds talk
 * to the Communicator (the singleton, to it's many instances).
 *
 * A User might use multiple clients and do different things on different machines
 *
 * A Client might be used by a different humans, and one human might be using
 * several Persons
 *
 * So, a Mind is where things get placed as being specific to one Person and one Client
 *
 * All communication between a client and a user goes through a Mind (it will
 * contextualize the session using its namespace, so this uuid is what goes into queries as a hint).
 * 
 * It also maintains the physical configuration of the client and persists that.
 *
 * @author John Christopher Anderson
 */
public class Mind 
		extends ObjBase {
	
    static final Logger logger = LoggerFactory.getLogger(Mind.class);
    
    /**
     * This is not a serverside API method
     * @param client
     * @param person
     * @return 
     */
    public static Mind create(Person person) {
        /* TODO check if this Mind already exists */
        //Create the new Mind object
        Mind out = new Mind();
        out.person = person;
        //out.save();
        return out;
    }

    public Mind() {
        engine = getEngine();
        userDictionary = new HashMap<>();
    }

    /**
     * Serverside API relayed method to 'clear the mind'
     */
    public final void clear() {
        engine = null;
        getEngine();
    }

    /* Similar to runCommand but does not "defuzzify".
     * Used for the "serverEval" channel. */
    public synchronized Object eval(String cmd, ScriptAPI api) throws ScriptException {
            return getEngine().eval(cmd, Language.JAVASCRIPT, api);    
    }

    /**
     * Try to run the command as if it were a proper script.  If fails,
     * try it again after injecting some objects that might not be in 
     * the engine yet.  If all that fails, return, then abort and disambiguate.
     * @param socket_id
     * @param cmd
     * @return 
     */
    /* TODO: race condition. ScriptEngine execution needs to be serialized. */
    public synchronized Object runCommand(String cmd, ScriptAPI api) throws ScriptException {
        return eval(cmd, api);
//                runCommandWithNamespace(cmd);
    } 

    /**
     * Return the scripting engine
     * (or a new instance if it doesn't already exist)
     */
    public MetaEngine getEngine() {
        if (engine == null) {
            engine = new MetaEngine();
        }
        return engine;
    }
    
    public Person getPerson() {
    	return this.person;
    }
    
    public List<Message> getLastCommands() {
        return this.lastCommands;
    }

    public void addLastCommand(Message command){
        lastCommands.add(command);
    }
    
    public void addLastCommand(Channel channel, Object data) {
        lastCommands.add(new Message(channel, data, null, null));
    }
    
    private Person person;
    
    @Setter
    private String username;
    
    private transient List<Message> lastCommands = new ArrayList<>();
    
    private MetaEngine engine;
    
    private  List<String> lastSharables = new ArrayList<>();
    @Getter @Setter 
    private transient ClientConnection connection;
    
    private Map<Date, String> commandHistory;
    
    private static final int MAX_SIZE_LAST_SHARABLES = 50;




    /* Maps "Doo ID" to Doo object */
    
    //JCA:  THIS IS PROBLEMATIC, IT IS HOLDING STATE FOR A DOO, NOT LOOSE-COUPLED
    //NOT SURE WHY THIS WOULD BE HERE INSTEAD OF THE HOPPER, I THINK THIS IS PROBABLY WRONG

    public ClientConnection getClientConnection() {
       return this.connection;
   }

    public String pullUUIDFromNamespace(String str) {
        System.out.println("JCA:  Need to implement mapping names/namespace tokens to UUIDs of Sharables");
        return null;
    }

    public List<String> getRecentSharables() {
        return this.lastSharables;
    }

    public void addLastSharable(String lastid) {
        if(this.lastSharables.contains(lastid)) {
            lastSharables.remove(lastid);
        }
        lastSharables.add(lastid);
        if(lastSharables.size() > MAX_SIZE_LAST_SHARABLES) {
            lastSharables.remove(MAX_SIZE_LAST_SHARABLES);
        }
    }
    
    //TODO: clean up interface here
    public Object evalFunction(String code, String name, List args, ScriptAPI api) throws ScriptException {
        return getEngine().invoke(code, name, args, api);
    }

    public Object invoke(Function f, List args, ScriptAPI api) throws ScriptException, NoSuchMethodException {
        return getEngine().invoke(f, convertArgs(args), api);
    }
    
    public Object invokeMethod(Module m, String methodName, List args, ScriptAPI api) throws NoSuchMethodException, ScriptException{
        return getEngine().invoke(m, methodName, convertArgs(args), api);
    }
    
    private static List convertArgs(List args){
        List out = new ArrayList();
        for (Object o : args){
            out.add(ScriptAPI.convertToNative(o));
        }
        return out;
    }

    private Map<String, ObjectId> userDictionary;
}
