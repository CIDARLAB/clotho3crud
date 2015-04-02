package org.clothocad.core.communication.ws;

import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.util.JSON;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;

import lombok.extern.slf4j.Slf4j;

import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.Subject;

import org.eclipse.jetty.websocket.WebSocket;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;

@Slf4j
public class ClothoWebSocket
        extends ClientConnection
        implements WebSocket.OnTextMessage {

    private WebSocket.Connection connection;
    private Subject subject;
    private final Router router;

    private class CallRouter implements Callable {

        private final Message message;
        private final ClientConnection connection;

        CallRouter(ClientConnection connection, Message message) {
            this.message = message;
            this.connection = connection;
        }

        @Override
        public Object call() throws Exception {
            router.receiveMessage(connection, message);
            return null;
        }
    }

    public ClothoWebSocket(String id, Router router) {
        super(id);
        this.router = router;
    }

    @Override
    public void onClose(int closeCode, String message) {
    }

    @Override
    public void send(Message msg) {
        try {
            String messageString = JSON.serializeForExternal(msg);
            connection.sendMessage(messageString);
            log.trace("Sent: {}", messageString);
        } catch (IOException ex) {
            log.error("Cannot send message", ex);
        }
    }

    @Override
    public void onMessage(String messageString) {
        log.trace("WebSocket #{} received message {}", this.getId(), messageString);
        try {
            Message message = JSON.mapper.readValue(messageString, Message.class);
            subject.execute(new CallRouter(this, message));
        } catch (JsonParseException ex) {
            log.error("WebSocket #{} received malformed message: {}", this.getId(), messageString);
        } catch (JsonMappingException ex) {
            throw new RuntimeException(ex);
        } catch (IOException ex) {
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }

    public boolean isOpen() {
        return connection.isOpen();
    }

    @Override
    public void onOpen(Connection connection) {
        log.debug("New connection opened. Connection id is {}.", this.getId());
        this.connection = connection;
        // Close WebSocket after 1 week ofidle time
        connection.setMaxIdleTime(7 * 24 * 3600000);

        subject = SecurityUtils.getSubject();

    }

    @Override
    public void deregister(Channel channel, String requestId) {
        String messageString = String.format("{\"channel\": \"%s\", \"requestId\": \"%s\"}", channel, requestId);
        try {
            connection.sendMessage(messageString);
            log.trace("Sent deregister: {}", messageString);
        } catch (IOException ex) {
            log.error("Cannot send deregister message:{}", ex.getMessage());
        }
    }
}
