/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public interface CredentialStore {
    
    public ClothoAccount getAccount(String username);
    public void saveAccount(ClothoAccount account);

    public AuthGroup getGroup(String groupName);        
    public void saveGroup(AuthGroup authGroup);

    public Map<String, Set<ClothoAction>> getUserPermissions(ObjectId id);
    public Map<String, Set<ClothoAction>> getGroupPermissions(ObjectId id);

    public void deleteAllCredentials();
}
