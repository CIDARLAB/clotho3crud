/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security.nosecurity;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.authz.Permission;
import org.apache.shiro.realm.Realm;
import org.apache.shiro.session.Session;
import org.apache.shiro.session.mgt.DefaultSessionManager;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.subject.Subject;
import org.apache.shiro.subject.SubjectContext;
import org.apache.shiro.web.mgt.DefaultWebSecurityManager;
import org.apache.shiro.web.mgt.DefaultWebSubjectFactory;
import org.apache.shiro.web.subject.WebSubjectContext;

/**
 *
 * @author spaige
 */
public class PermissiveSecurityManager extends DefaultWebSecurityManager {

    public PermissiveSecurityManager() {
        super();
        setSubjectFactory(new AuthoringSubjectFactory());
        setSessionManager(new DefaultSessionManager());
    }

    @SuppressWarnings({"UnusedDeclaration"})
    public PermissiveSecurityManager(Realm singleRealm) {
        this();
        setRealm(singleRealm);
    }

    @SuppressWarnings({"UnusedDeclaration"})
    public PermissiveSecurityManager(Collection<Realm> realms) {
        this();
        setRealms(realms);
    }

    @Override
    public void checkPermission(PrincipalCollection principals, Permission permission) throws AuthorizationException {
    }

    @Override
    public void checkPermission(PrincipalCollection principals, String permission) throws AuthorizationException {
    }

    @Override
    public void checkPermissions(PrincipalCollection principals, Collection<Permission> permissions) throws AuthorizationException {
    }

    @Override
    public void checkPermissions(PrincipalCollection principals, String... permissions) throws AuthorizationException {
    }

    @Override
    public boolean[] isPermitted(PrincipalCollection principals, List<Permission> permissions) {
        boolean[] trues = new boolean[permissions.size()];
        Arrays.fill(trues,true);
        return trues;
    }

    @Override
    public boolean isPermitted(PrincipalCollection principals, Permission permission) {
        return true;
    }

    @Override
    public boolean isPermitted(PrincipalCollection principals, String permissionString) {
        return true;
    }

    @Override
    public boolean[] isPermitted(PrincipalCollection principals, String... permissions) {
        boolean[] trues = new boolean[permissions.length];
        Arrays.fill(trues,true);
        return trues;       
    }

    @Override
    public boolean isPermittedAll(PrincipalCollection principals, Collection<Permission> permissions) {
        return true;
    }

    @Override
    public boolean isPermittedAll(PrincipalCollection principals, String... permissions) {
        return true;
    }
    public static class AuthoringSubjectFactory extends DefaultWebSubjectFactory {

        @Override
        public Subject createSubject(SubjectContext context) {
            if (!(context instanceof WebSubjectContext)) {
                org.apache.shiro.mgt.SecurityManager securityManager = context.resolveSecurityManager();
                Session session = context.resolveSession();
                boolean sessionCreationEnabled = context.isSessionCreationEnabled();
                PrincipalCollection principals = context.resolvePrincipals();
                boolean authenticated = context.resolveAuthenticated();
                String host = context.resolveHost();

                return new LoggedInSubject(principals, authenticated, host, session, sessionCreationEnabled, securityManager);
            }
            WebSubjectContext wsc = (WebSubjectContext) context;
            org.apache.shiro.mgt.SecurityManager securityManager = wsc.resolveSecurityManager();
            Session session = wsc.resolveSession();
            boolean sessionEnabled = wsc.isSessionCreationEnabled();
            PrincipalCollection principals = wsc.resolvePrincipals();
            boolean authenticated = wsc.resolveAuthenticated();
            String host = wsc.resolveHost();
            ServletRequest request = wsc.resolveServletRequest();
            ServletResponse response = wsc.resolveServletResponse();

            return new LoggedInWebSubject(principals, authenticated, host, session, sessionEnabled,
                    request, response, securityManager);
        }
    }
}