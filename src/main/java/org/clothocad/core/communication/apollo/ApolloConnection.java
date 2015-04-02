package org.clothocad.core.communication.apollo;

import java.util.UUID;

import javax.jms.Session;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.apollo.ClothoMessageProducer;

import org.clothocad.core.communication.ClientConnection;

public class ApolloConnection 
	extends ClientConnection {

	private Session session;
	private String correlationId;
        private ClothoMessageProducer messageProducer;
	
	public ApolloConnection(String id, Session session, String correlationId) {
		super(id);
		this.session = session;
		this.correlationId = correlationId;
        this.messageProducer = new ClothoMessageProducer(this);
	}
	
	public ApolloConnection(Session session) {
		super(UUID.randomUUID().toString());
		this.session = session;
	}

	public Session getSession() {
		return this.session;
	}
	
	public String getCorrelationId() {
		return this.correlationId;
	}

    @Override
    public void send(Message msg) {
        messageProducer.onSuccess(msg.getData());
    }
}
