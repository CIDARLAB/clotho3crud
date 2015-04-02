package org.clothocad.core.communication.apollo;

import java.io.IOException;
import java.util.Map;
import java.util.UUID;

import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.Session;
import org.clothocad.core.util.JSON;
import org.fusesource.stomp.jms.message.StompJmsMessage;

public class ClothoMessageConsumer
        implements Runnable {

    private Session session;
    private StompJmsMessage message;

    public ClothoMessageConsumer(Session session, Message message)
            throws Exception {
        this.session = session;

        if (message != null && message instanceof StompJmsMessage) {
            this.message = (StompJmsMessage) message;
        } else {
            throw new Exception("INVALID MESSAGE!");
        }
    }

    @Override
    public void run() {
        try {
            if (this.message.propertyExists("request")) {

                Map<String, Object> json = JSON.deserializeObjectToMap(message.getStringProperty("request"));

                // get the message's correlation id
                String sCorrelationID = this.message.getJMSCorrelationID();
                if (null == sCorrelationID) {
                    sCorrelationID = UUID.randomUUID().toString();
                    this.message.setJMSCorrelationID(sCorrelationID);
                }


                String auth_key = (String) json.get("auth_key");
                if (auth_key != null) {
                    //XXX:
                    //auth_key will not be stored in the body of the request from the clotho client
                    //an equivalent in the Stomp layer is fine
                    //whatever value JMS uses to uniquely identify endpoints, we should use that to retrieve the appropriate connection and submit that to the router along with the message
                    /*ApolloConnection connection = null;
                    Mind mind = Router.get().getMind(auth_key);
                    if (null == mind.getClientConnection()) {
                        connection = new ApolloConnection(
                                session);
                    } else {
                        ClientConnection conn = mind.getClientConnection();
                        if (conn instanceof ApolloConnection) {
                            connection = (ApolloConnection) conn;
                        } else {
                            Router.get().receiveMessage(
                                    connection,
                                    json.getString("channel"),
                                    json);
                            return;
                        }
                    }

                    // store the callback-handler in the callback-handler table
                    CallbackHandlerTable.put(
                            connection.getId(),
                            new CallbackHandler(this.session, this.message));

                    // route the message
                    Router.get().receiveMessage(
                            connection, new org.clothocad.core.layers.communication.Message(Channel.valueOf(json.get("channel").toString()), json));*/
                } else {
                    // this needs to be clarified...
                }
            }
        } catch (JMSException e) {
            e.printStackTrace();
        }  catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}
