/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.apache.shiro.authz.AuthorizationInfo;
import org.apache.shiro.authz.Permission;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public class ClothoAuthorizationInfo extends PermissionHolder implements AuthorizationInfo {
    
    private Set<String> groups;

    public ClothoAuthorizationInfo(Set<String> groups, Map<ObjectId, PermissionsOnObject> objectPermissions) {
        this.groups = groups;
        this.permissions = objectPermissions;
    }

    public ClothoAuthorizationInfo(){
        groups = new HashSet<>();
        permissions = new HashMap<>();
    }
    
    @Override
    public Collection<String> getRoles() {
        return groups;
    }

    @Override
    public Collection<String> getStringPermissions() {
        Collection<String> out = new HashSet<String>();
        for (PermissionsOnObject p : permissions.values()){
            out.addAll(p.permissions);
        }
        return out;
    }

    @Override
    public Collection<Permission> getObjectPermissions() {
        return Collections.emptySet();
    }
    
    public Collection<String> getStringPermissions(ObjectId id){
        return permissions.get(id).permissions;
    }
    
    public void addGroup(String group){
        groups.add(group);
    }
    
    public void removeGroup(String group){
        groups.remove(group);
    }

}
