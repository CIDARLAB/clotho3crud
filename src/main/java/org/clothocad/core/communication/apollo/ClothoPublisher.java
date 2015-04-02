package org.clothocad.core.communication.apollo;

import java.util.Map;
import javax.jms.Connection;
import javax.jms.DeliveryMode;
import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.fusesource.stomp.jms.StompJmsConnectionFactory;
import org.fusesource.stomp.jms.StompJmsDestination;

public class ClothoPublisher {

	private Connection connection = null;
	private Session session = null;
	private MessageProducer producer = null;
	
	public ClothoPublisher() 
			throws Exception {
		StompJmsConnectionFactory factory = 
				new StompJmsConnectionFactory();
		factory.setBrokerURI("tcp://localhost:61613");
		
		this.connection = factory.createConnection("admin", "password");
		this.connection.start();
		
		// third, we create a session for the connection
        this.session = connection.createSession(
        		false, Session.CLIENT_ACKNOWLEDGE);
        
        this.producer = session.createProducer(
        		new StompJmsDestination("/topic/CLOTHO.UPDATE"));
        this.producer.setDeliveryMode(DeliveryMode.PERSISTENT);

	}
	
	public void publish(Map<String, Object> json) {
		try {
			System.out.println("[ClothoPublisher.publish] -> "+json);
			
			Message response = session.createMessage();
			response.setStringProperty("RESPONSE", json.toString()); //XXX: this class is pretty out of date
			this.producer.send(response);
		} catch (JMSException e) {
			e.printStackTrace();
		}
	}
}
