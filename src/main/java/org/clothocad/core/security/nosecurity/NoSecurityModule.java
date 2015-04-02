/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security.nosecurity;

import org.clothocad.core.security.nosecurity.PermissiveSecurityManager;
import com.google.inject.binder.AnnotatedBindingBuilder;
import java.util.Collection;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.config.ConfigurationException;
import org.apache.shiro.web.mgt.WebSecurityManager;
import org.clothocad.core.security.SecurityModule;
import org.eclipse.jetty.servlet.ServletContextHandler;

/**
 *
 * @author spaige
 */
public class NoSecurityModule extends SecurityModule {

    public NoSecurityModule() {
        super();
    }

    public NoSecurityModule(ServletContextHandler servletContextHandler) {
        super(servletContextHandler);
    }

    
    @Override
    protected void bindWebSecurityManager(AnnotatedBindingBuilder<? super WebSecurityManager> bind) {
        try {
            // bind.toConstructor(PermissiveSecurityManager.class.getConstructor(Collection.class)).asEagerSingleton();
            bind.toConstructor(PermissiveSecurityManager.class.getConstructor(Collection.class)).asEagerSingleton();
        } catch (NoSuchMethodException e) {
            throw new ConfigurationException("This really shouldn't happen.  Either something has changed in Shiro, or there's a bug in ShiroModule.", e);
        }
    }
}
