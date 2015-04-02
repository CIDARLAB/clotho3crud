package org.clothocad.core.communication.apollo;

import java.util.HashMap;
import java.util.Map;
import org.clothocad.core.communication.Callback;

public class ClothoMessageProducer 
		implements Callback {
	
	private ApolloConnection connection;
	public ClothoMessageProducer(ApolloConnection connection) {
		this.connection = connection;
	}

	@Override
	public void onSuccess(Object json) {
            CallbackHandler cbh = 
                    CallbackHandlerTable.get(connection.getId());
            if(null != cbh) {
                System.err.println("[ClothoMessageProducer.onSuccess] -> "+json+" -> "+connection.getId());
                cbh.respond(json);			
            }
	}

	@Override
	public void onFailure(Throwable err) {
            CallbackHandler cbh = CallbackHandlerTable.get(
                            connection.getId());
            if(null != cbh) {
                Map<String,Object> jsonResponse = new HashMap<>();
                try {
                    jsonResponse.put("CLOTHO-ERROR", err.getLocalizedMessage());
                } catch(Exception e) {}
                System.err.println("[ClothoMessageProducer.onFailure] -> "+jsonResponse+" -> "+connection.getId());
                cbh.respond(jsonResponse);			
            }
	}

}
