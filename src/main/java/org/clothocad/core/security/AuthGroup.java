/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.HashMap;
import java.util.Map;
import lombok.Getter;
import org.clothocad.core.datums.ObjectId;
import org.jongo.marshall.jackson.oid.Id;

/**
 *
 * @author spaige
 */
public class AuthGroup extends PermissionHolder {
    
    @JsonCreator
    public AuthGroup(@JsonProperty("name") String name, @JsonProperty("permissions") Map<ObjectId, PermissionsOnObject> permissions){
        this.name = name;
        this.permissions = permissions;
    }
    
    public AuthGroup(String name){
        this.name = name;
        this.permissions = new HashMap<>();
    }
    
    @Getter
    @Id
    protected String name;

}
