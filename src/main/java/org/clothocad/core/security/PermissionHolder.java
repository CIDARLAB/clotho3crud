/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import lombok.Getter;
import org.apache.shiro.util.CollectionUtils;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public abstract class PermissionHolder {
    
    @Getter
    protected Map<ObjectId, PermissionsOnObject> permissions;
    
    //String permission == Shiro wildcard permssion string
    public static ObjectId getObjectId(String permission) {
        List<String> parts = CollectionUtils.asList(permission.split(":"));
        return new ObjectId(parts.get(2));
    }
    
    public static ClothoAction getAction(String stringPermission){
        return ClothoAction.valueOf(CollectionUtils.asList(stringPermission.split(":")).get(1));
    }

    public static String constructPermissionString(ClothoAction permission, ObjectId id){
        return "data:" + permission.name() + ":" + id.toString();
    }
    
    public void addPermission(ClothoAction permission, ObjectId id){
        if (!permissions.containsKey(id)){
            permissions.put(id, new PermissionsOnObject());
        }
        permissions.get(id).permissions.add(constructPermissionString(permission, id));
    }
    
    public void removePermission(ClothoAction permission, ObjectId id){
        if (permissions.containsKey(id)){
            permissions.get(id).permissions.remove(constructPermissionString(permission, id));
        }
    }    
    
    public Set<ClothoAction> getActions(ObjectId id){
        Set<ClothoAction> actions = EnumSet.noneOf(ClothoAction.class);
        
        for(String stringPermission : permissions.get(id).permissions){
            actions.add(getAction(stringPermission));
        }
        
        return actions;
    }
    
}
