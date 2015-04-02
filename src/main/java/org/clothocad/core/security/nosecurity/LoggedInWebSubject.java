/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security.nosecurity;

import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.mgt.SecurityManager;
import org.apache.shiro.session.Session;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.web.subject.support.WebDelegatingSubject;
import static org.clothocad.core.security.ServerSubject.SERVER_USER;

/**
 *
 * @author spaige
 */
public class LoggedInWebSubject extends WebDelegatingSubject {

    public LoggedInWebSubject(PrincipalCollection principals, boolean authenticated, String host, Session session, ServletRequest request, ServletResponse response, SecurityManager securityManager) {
        super(principals, authenticated, host, session, request, response, securityManager);
    }

    public LoggedInWebSubject(PrincipalCollection principals, boolean authenticated, String host, Session session, boolean sessionEnabled, ServletRequest request, ServletResponse response, SecurityManager securityManager) {
        super(principals, authenticated, host, session, sessionEnabled, request, response, securityManager);
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
