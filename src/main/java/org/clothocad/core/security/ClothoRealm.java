/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import javax.inject.Inject;
import javax.inject.Singleton;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.authc.AccountException;
import org.apache.shiro.authc.AuthenticationException;
import org.apache.shiro.authc.AuthenticationInfo;
import org.apache.shiro.authc.AuthenticationToken;
import org.apache.shiro.authc.SimpleAuthenticationInfo;
import org.apache.shiro.authc.UsernamePasswordToken;
import org.apache.shiro.authc.credential.HashedCredentialsMatcher;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.authz.AuthorizationInfo;
import org.apache.shiro.authz.permission.RolePermissionResolver;
import org.apache.shiro.crypto.hash.Sha256Hash;
import org.apache.shiro.realm.AuthorizingRealm;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.subject.SimplePrincipalCollection;
import org.apache.shiro.subject.Subject;
import org.clothocad.core.datums.ObjectId;
import static org.clothocad.core.security.ServerSubject.SERVER_USER;

/**
 *
 * Note: probably worth implementing a more intelligent permission search for
 * data object permissions
 * 
 * TODO: locking
 * TODO: credential store throws entitynotfound exception
 *
        //TODO: Networking: realm name is clotho instance name
        //                  add new user to local instance group 
 * @author spaige
 */
@Slf4j
@Singleton
public class ClothoRealm extends AuthorizingRealm {
    
    public static final String ALL_GROUP = "_all";
    public static final String ANONYMOUS_USER = "_anonymous";
    
    private CredentialStore store;
    
    public static final AuthenticationToken getAnonymousUserToken(){
        return new UsernamePasswordToken(ANONYMOUS_USER, ANONYMOUS_USER);
    }

    @Inject
    public ClothoRealm(CredentialStore store, RolePermissionResolver roleResolver) {
        super();
        
        //XXX: up number of iterations
        HashedCredentialsMatcher matcher = new HashedCredentialsMatcher(Sha256Hash.ALGORITHM_NAME);
        matcher.setStoredCredentialsHexEncoded(false);
        
        this.store = store;
        setAuthenticationTokenClass(UsernamePasswordToken.class);
        setCredentialsMatcher(matcher);

        setRolePermissionResolver(roleResolver);

        setUpRealm();
    }
    
    protected void setUpRealm(){
        //set up default groups
        if (store.getGroup(ALL_GROUP) == null){
            addGroup(ALL_GROUP);
        }

        //set up anonymous user
        if (store.getAccount(ANONYMOUS_USER) == null){
            addAccount(ANONYMOUS_USER, ANONYMOUS_USER);
        }
        
        if (store.getAccount(SERVER_USER) == null){
            store.saveAccount(new DummyAccount(SERVER_USER));
        }
    }

    @Override
    protected Object getAuthorizationCacheKey(PrincipalCollection principals) {
        return principals.fromRealm(getRealmName()).iterator().next();
    }
    
    @Override
    protected AuthorizationInfo doGetAuthorizationInfo(PrincipalCollection pc) {
        log.debug("getting authz info for {}", pc);

        ClothoAccount account = store.getAccount(pc.getPrimaryPrincipal().toString());
        return account.getAuthzInfo();
    }

    @Override
    protected AuthenticationInfo doGetAuthenticationInfo(AuthenticationToken at) throws AuthenticationException {
        log.debug("getting authc info for {}", at);
        
        ClothoAccount account = store.getAccount(((UsernamePasswordToken) at).getUsername());
        if (!account.isAuthenticatable()) throw new AccountException("Cannot authenticate as " + at.getPrincipal().toString());
        return account;
    }
    
    public void addAccount(String username, String password) {
        if (store.getAccount(username) != null){
            throw new javax.persistence.EntityExistsException();
        }
        store.saveAccount(new ClothoAccount(username, password));
    }
        
        
    public void addGroup(String groupName){
        if (store.getAccount(groupName) != null){
            throw new javax.persistence.EntityExistsException();
        }
        store.saveGroup(new AuthGroup(groupName));
    }
    

    public void addPermission(String username, ClothoAction permission, ObjectId id) {
        addPermission(username, permission, id, true);
    }
    
    public void addPermission(String username, ClothoAction permission, ObjectId id, boolean checkPriv) {
        if (checkPriv) checkCurrentSubjectGrant(id);
        ClothoAccount account = store.getAccount(username);
        account.getAuthzInfo().addPermission(permission, id);
        store.saveAccount(account);
        if (isAuthorizationCachingEnabled()) clearCachedAuthorizationInfo(new SimplePrincipalCollection(username, getRealmName()));
    }

    public void addPermissionToGroup(String groupName, ClothoAction permission, ObjectId id) {
        checkCurrentSubjectGrant(id);
        AuthGroup group = store.getGroup(groupName);
        group.addPermission(permission, id);
        store.saveGroup(group);
        if (getAuthorizationCache() != null) clearCachedAuthorizationInfo(new SimplePrincipalCollection(groupName, getRealmName()));
    }

    public void removePermissionFromGroup(String groupName, ClothoAction permission, ObjectId id) {
        checkCurrentSubjectGrant(id);
        AuthGroup group = store.getGroup(groupName);
        group.removePermission(permission, id);
        store.saveGroup(group);
        if (getAuthorizationCache() != null) clearCachedAuthorizationInfo(new SimplePrincipalCollection(groupName, getRealmName()));
    }

    public void removePermission(String username, ClothoAction permission, ObjectId id) {        
        checkCurrentSubjectGrant(id);
        ClothoAccount account = store.getAccount(username);
        account.getAuthzInfo().removePermission(permission, id);
        store.saveAccount(account);
        if (getAuthorizationCache() != null) clearCachedAuthorizationInfo(new SimplePrincipalCollection(username, getRealmName()));
    }

    public void removePermissions(String username, Set<ClothoAction> permissions, ObjectId id){
        for (ClothoAction permission: permissions) {
            removePermission(username, permission, id);
        }
    }
    
    public void addPermissions(String username, Set<ClothoAction> permissions, Set<ObjectId> ids) {
        for (ClothoAction permission : permissions) {
            for (ObjectId id : ids) {
                try {
                    addPermission(username, permission, id);
                } catch (AuthorizationException e) {
                    //log 
                }
            }
        }
    }
    public void addPermissions(String username, Set<ClothoAction> permissions, ObjectId id) {
        addPermissions(username, permissions, id, true);
    }
    
    
    public void addPermissions(String username, Set<ClothoAction> permissions, ObjectId id, boolean checkPriv) {
        for (ClothoAction permission : permissions) {
            addPermission(username, permission, id, checkPriv);
        }
    }

    public void addPermissionsToGroup(String groupName, Set<ClothoAction> permissions, ObjectId id){
        for (ClothoAction permission : permissions){
            addPermissionToGroup(groupName, permission, id);
        }
    }
    
    public void removePermissionsFromGroup(String groupName, Set<ClothoAction> permissions, ObjectId id){
        for (ClothoAction permission : permissions){
            removePermissionFromGroup(groupName, permission, id);
        }
    }
    
    public void deleteAll() {
        store.deleteAllCredentials();
        setUpRealm();
    }

    public void removePublic(ObjectId id) {
        checkCurrentSubjectGrant(id);
        removePermissionsFromGroup(ALL_GROUP, ClothoPermission.READ.actions ,id);
    }

    public void setPublic(ObjectId id) {
        checkCurrentSubjectGrant(id);
        addPermissionsToGroup(ALL_GROUP, ClothoPermission.READ.actions, id);
    }
    
    public String getRealmName() {
        return "clotho";
    }
    
    private static void checkCurrentSubjectGrant(ObjectId id) throws AuthorizationException {
        SecurityUtils.getSubject().checkPermission("data:grant:"+id.toString());
    }
    
    private static void checkCurrentSubjectIs(String username) throws AuthorizationException {
        if (!SecurityUtils.getSubject().getPrincipal().equals(username))
            throw new AuthorizationException("Current subject is not "+ username);
    }
    
    public void updatePassword(String username, String password) { 
        ClothoAccount account = store.getAccount(username);
        account.setPassword(password);
        store.saveAccount(account);
    }
    
    public void addPrincipal(String username, Object principal, String realm){
        checkCurrentSubjectIs(username);
        ClothoAccount account = store.getAccount(username);
        account.merge(new SimpleAuthenticationInfo(principal, null, realm));
        store.saveAccount(account);
    }
    
    public Map<String, Set<ClothoAction>> getUserPermissions(ObjectId id){
        checkCurrentSubjectGrant(id);
        return store.getUserPermissions(id);
    }
    
    public Map<String, Set<ClothoAction>> getGroupPermissions(ObjectId id){
        checkCurrentSubjectGrant(id);
        return store.getGroupPermissions(id);
    }
    
    public Set<ClothoAction> getCurrentSubjectPermissions(ObjectId id){
        Subject currentSubject = SecurityUtils.getSubject();
        currentSubject.checkPermission("data:view:"+id.toString());
        
        Set<ClothoAction> permittedActions = new HashSet<>();
        for (ClothoAction action : ClothoAction.values()){
            if (currentSubject.isPermitted("data:"+action.name()+":"+id.toString()))
                permittedActions.add(action);
        }
        
        return permittedActions;
    }
}
