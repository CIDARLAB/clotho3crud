/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.google.inject.AbstractModule;
import org.apache.shiro.authz.permission.RolePermissionResolver;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.security.CredentialStore;

/**
 *
 * @author spaige
 */
public class JongoModule extends AbstractModule {

    @Override
    protected void configure() {
        bind(ClothoConnection.class).to(JongoConnection.class);
        bind(CredentialStore.class).to(JongoConnection.class);
        bind(RolePermissionResolver.class).to(JongoConnection.class);
    }
    
}
