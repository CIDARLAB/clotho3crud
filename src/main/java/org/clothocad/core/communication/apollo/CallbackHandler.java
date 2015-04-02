package org.clothocad.core.communication.apollo;

import javax.jms.DeliveryMode;
import javax.jms.Message;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.fusesource.stomp.jms.StompJmsDestination;

public class CallbackHandler {

	private Session session;
	private Message request;
	private MessageProducer replyProducer;
	
	public CallbackHandler(Session session, Message request) {
		this.session = session;
		this.request = request;
		
		try {
			// Setup a message producer to respond to messages from clients, we will get the destination
	        // to send to from the JMSReplyTo header field from a Message
	        this.replyProducer = this.session.createProducer(
	        		new StompJmsDestination(ClothoConstants.CLOTHO_RESPONSE_QUEUE));
	        this.replyProducer.setDeliveryMode(DeliveryMode.PERSISTENT);

		} catch(Exception e) {
			e.printStackTrace();
		}
    }
	
        //XXX: this probably needs to be revisited, use JSON.serialize, etc
	public void respond(Object json) {
		try {
			Message response = session.createMessage();
			response.setJMSCorrelationID(request.getJMSCorrelationID());
			response.setStringProperty("response", json.toString());

                        System.out.println("[CallbackHandler.response] -> "+request.getJMSCorrelationID()
                                +" -> "+ClothoConstants.CLOTHO_RESPONSE_QUEUE);
			this.replyProducer.send(response);
		} catch(javax.jms.InvalidDestinationException e) {
			// something went wrong while responding to the client...
			e.printStackTrace();
		} catch(javax.jms.JMSException e) {
			// something went wrong with the message ...
			e.printStackTrace();
		}
	}
}
