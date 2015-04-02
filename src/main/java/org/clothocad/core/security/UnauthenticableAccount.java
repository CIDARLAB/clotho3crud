/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import org.apache.shiro.authc.SimpleAuthenticationInfo;

/**
 * For system accounts that should not be logged in as by clients
 * @author spaige
 */
public class UnauthenticableAccount extends ClothoAccount {
    public UnauthenticableAccount(String username){
        super(username);
    }

    @JsonCreator
    public UnauthenticableAccount(@JsonProperty("_id") String username,
    @JsonProperty("authzInfo") ClothoAuthorizationInfo authzInfo) {
        super(username, null, authzInfo);
    
    }
    
    @Override
    public boolean isAuthenticatable() {
        return false;
    }
    
    
}
