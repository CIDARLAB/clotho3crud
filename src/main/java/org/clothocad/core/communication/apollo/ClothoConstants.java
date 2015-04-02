package org.clothocad.core.communication.apollo;

public class ClothoConstants {
	public static final String SERVER_URL = "tcp://localhost:61613";	
        public static final String CLOTHO_QUEUE = "/queue/CLOTHO";
        public static final String CLOTHO_RESPONSE_QUEUE = "/queue/CLOTHORESPONSE";
        
	public static final String CHANNEL = "channel";
	public static final String ACTION = "action";
	public static final String DATA = "data";
	public static final String CORRELATION_ID = "correlation";
	public static final String AUTHENTICATION = "authentication";
        
        // here we could also add some more constants, for example, indicating the JSON field names
}
