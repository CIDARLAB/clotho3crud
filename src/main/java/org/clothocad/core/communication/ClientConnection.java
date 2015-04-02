package org.clothocad.core.communication;

public abstract class ClientConnection {

    private String id;

    public ClientConnection(String id) {
        this.id = id;
    }

    public String getId() {
        return this.id;
    }

    public abstract void send(Message msg);

    public void deregister(Channel channel, String requestId) {
    }
}
