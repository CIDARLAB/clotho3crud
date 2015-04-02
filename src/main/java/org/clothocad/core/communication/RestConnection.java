/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import java.util.ArrayList;
import java.util.List;
import java.io.IOException;

import org.clothocad.core.util.JSON;

/**
 *
 * @author mosnicholas
 */
public class RestConnection extends ClientConnection {

    private boolean done = false;
    private String jsonResult;

    public RestConnection(String id) {
        super(id);
    }
    
    @Override
    public void send(Message msg) {
        try {
            jsonResult = JSON.serialize(msg);
            this.done = true;
        } catch (IOException e) {
            jsonResult = e.toString();
            // jsonResult = "Error";
        }
    }
    
    public boolean isDone() {
        return this.done;
    }

    public String getResult() {
        return jsonResult;
    }
}
