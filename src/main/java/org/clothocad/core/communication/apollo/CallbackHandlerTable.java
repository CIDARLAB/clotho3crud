package org.clothocad.core.communication.apollo;

import java.util.Hashtable;

public class CallbackHandlerTable {

	// key   ... the message's correlation id
	// value ... the callback-handler object 
	private static Hashtable<String, CallbackHandler> htCallbackHandlers;
	
	public static void put(String sConnectionID, CallbackHandler cbh) {
            if(null == htCallbackHandlers) {
                htCallbackHandlers = new Hashtable<String, CallbackHandler>();
            }
		
            if(!htCallbackHandlers.containsKey(sConnectionID)) {
		htCallbackHandlers.put(sConnectionID, cbh);
            } else {
            }        
	}
	
	public static CallbackHandler get(String sConnectionID) {
            if(null != htCallbackHandlers) {
                CallbackHandler cbh = htCallbackHandlers.get(sConnectionID);
                htCallbackHandlers.remove(sConnectionID);
                return cbh;
            }
            return (CallbackHandler)null;
	}
}
