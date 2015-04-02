/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security.nosecurity;

import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.mgt.SecurityManager;
import org.apache.shiro.session.Session;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.subject.support.DelegatingSubject;
import static org.clothocad.core.security.ServerSubject.SERVER_USER;

/**
 *
 * @author spaige
 */
public class LoggedInSubject extends DelegatingSubject {

    public LoggedInSubject(PrincipalCollection principals, boolean authenticated, String host, Session session, boolean sessionCreationEnabled, SecurityManager securityManager) {
        super(principals, authenticated, host, session, sessionCreationEnabled, securityManager);
    }
    
    @Override
    public Object getPrincipal() {
        return SERVER_USER;
    }
    
    @Override
    public boolean isAuthenticated() {
        return true;
    }

    @Override
    protected void assertAuthzCheckPossible() throws AuthorizationException {
    }
}
