/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import lombok.Getter;
/**
 * Please note this class uses Shiro wildcard-formatted string permissions, instead of 
 * ClothoPermission/ClothoAction like other classes
 * 
 * @author spaige
 */
public class PermissionsOnObject {

    public PermissionsOnObject() {
        permissions = new HashSet<>();
        fieldPermissions = new HashMap<>();
    }

    @JsonCreator
    public PermissionsOnObject(@JsonProperty("permissions") Set<String> permissions, 
    @JsonProperty("fieldPermissions") Map<String, Set<String>> fieldPermissions) {
        this.permissions = permissions;
        this.fieldPermissions = fieldPermissions;
    }

    @Getter
    protected Set<String> permissions;
    
    
    protected Map<String, Set<String>> fieldPermissions;

}
