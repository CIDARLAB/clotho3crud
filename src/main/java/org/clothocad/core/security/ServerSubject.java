/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import org.apache.shiro.authc.AuthenticationException;
import org.apache.shiro.authc.AuthenticationToken;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.authz.Permission;
import org.apache.shiro.session.Session;
import org.apache.shiro.subject.ExecutionException;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.subject.Subject;
import org.apache.shiro.subject.support.SubjectCallable;
import org.apache.shiro.subject.support.SubjectRunnable;

/**
 *
 * @author spaige
 */
public class ServerSubject implements Subject {

    public static String SERVER_USER = "_server";
    
    private boolean isRunAs;
    private Object runningAs;
    
    @Override
    public Object getPrincipal() {
        if (isRunAs) return runningAs;
        return SERVER_USER;
    }

    @Override
    public PrincipalCollection getPrincipals() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean isPermitted(String permission) {
        return true;
    }

    @Override
    public boolean isPermitted(Permission permission) {
        return true;
    }

    @Override
    public boolean[] isPermitted(String... permissions) {
        boolean[] out = new boolean[permissions.length];
        for (int i = 0; i < permissions.length; i++){
            out[i] = true;
        }
        return out;
    }

    @Override
    public boolean[] isPermitted(List<Permission> permissions) {
        boolean[] out = new boolean[permissions.size()];
        for (int i = 0; i < permissions.size(); i++){
            out[i] = true;
        }
        return out;
    }

    @Override
    public boolean isPermittedAll(String... permissions) {
        return true;
    }

    @Override
    public boolean isPermittedAll(Collection<Permission> permissions) {
        return true;
    }

    @Override
    public void checkPermission(String permission) throws AuthorizationException {
    }

    @Override
    public void checkPermission(Permission permission) throws AuthorizationException {
    }

    @Override
    public void checkPermissions(String... permissions) throws AuthorizationException {
    }

    @Override
    public void checkPermissions(Collection<Permission> permissions) throws AuthorizationException {
    }

    @Override
    public boolean hasRole(String roleIdentifier) {
        return true;
    }

    @Override
    public boolean[] hasRoles(List<String> roleIdentifiers) {
        boolean[] out = new boolean[roleIdentifiers.size()];
        for (int i = 0; i < roleIdentifiers.size(); i++){
            out[i] = true;
        }
        return out;
    }

    @Override
    public boolean hasAllRoles(Collection<String> roleIdentifiers) {
        return true;
    }

    @Override
    public void checkRole(String roleIdentifier) throws AuthorizationException {
    }

    @Override
    public void checkRoles(Collection<String> roleIdentifiers) throws AuthorizationException {
    }

    @Override
    public void checkRoles(String... roleIdentifiers) throws AuthorizationException {
    }

    @Override
    public void login(AuthenticationToken token) throws AuthenticationException {
    }

    @Override
    public boolean isAuthenticated() {
        return true;
    }

    @Override
    public boolean isRemembered() {
        return false;
    }

    @Override
    public Session getSession() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Session getSession(boolean create) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void logout() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public <V> V execute(Callable<V> callable) throws ExecutionException {
        Callable<V> associated = associateWith(callable);
        try {
            return associated.call();
        } catch (Throwable t) {
            throw new ExecutionException(t);
        }
    }

    @Override
    public void execute(Runnable runnable) {
        Runnable associated = associateWith(runnable);
        associated.run();
    }

    @Override
    public <V> Callable<V> associateWith(Callable<V> callable) {
        return new SubjectCallable<>(this, callable);
    }

    @Override
    public Runnable associateWith(Runnable runnable) {
        if (runnable instanceof Thread) {
            String msg = "This implementation does not support Thread arguments because of JDK ThreadLocal " +
                    "inheritance mechanisms required by Shiro.  Instead, the method argument should be a non-Thread " +
                    "Runnable and the return value from this method can then be given to an ExecutorService or " +
                    "another Thread.";
            throw new UnsupportedOperationException(msg);
        }
        return new SubjectRunnable(this, runnable);
    }

    @Override
    public void runAs(PrincipalCollection principals) throws NullPointerException, IllegalStateException {
        isRunAs = true;
        runningAs = principals.getPrimaryPrincipal();
    }

    @Override
    public boolean isRunAs() {
        return isRunAs;
    }

    @Override
    public PrincipalCollection getPreviousPrincipals() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public PrincipalCollection releaseRunAs() {
        isRunAs = false;
        runningAs = null;
        return null;
    }
    
}
