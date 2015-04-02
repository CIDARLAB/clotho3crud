/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public class DummyAccount extends UnauthenticableAccount {

    
    @JsonCreator
    public DummyAccount( @JsonProperty("_id") String username) {
        super(username);
    }

    @Override
    public void addGroup(String group) {
    }

    @Override
    public void addPermission(ClothoAction permission, ObjectId id) {
    }

    @Override
    public void removeGroup(String group) {
    }

    @Override
    public void removePermission(ClothoAction permission, ObjectId id) {
    }

}
