package org.clothocad.core;

import com.google.common.collect.Lists;
import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Key;
import com.google.inject.Module;
import java.util.List;
import java.util.Properties;
import org.apache.shiro.SecurityUtils;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.nosecurity.NoSecurityModule;
import org.clothocad.model.BasicPart;
import org.clothocad.model.Part;
import org.apache.shiro.mgt.SecurityManager;
import org.clothocad.model.FreeForm;

/**
 *
 * @author spaige
 */
public class PersistorModularityDemo {
    public static void main(String[] args){

        
        //Setup: make a builder configured with the modules you want, then get 
        //the classes you want to use
        ClothoBuilder builder = new ClothoBuilder(
                //We have to have a security manager of some kind - this module 
                //sets up one that just does 'pass-through' checks on permissions
                new NoSecurityModule(),
                //Sets up the persistor to use Jongo/MongoDB as the database
                new JongoModule());
        
        //get a persistor
        Persistor p = builder.get(Persistor.class);
        
        //use the persistor to do stuff
        Part part = new BasicPart("Test Part", "this is a test part", "ATCG", new FreeForm(), null);
        
        p.save(part);
        
        Part retrievedPart = p.get(Part.class, part.getId());
    }
    
    public static class ClothoBuilder{
        ClothoBuilder(Properties config, Module... modulesArray){
            List<Module> modules = Lists.asList(new ClothoModule(config), modulesArray);
            injector = Guice.createInjector(modules);
            SecurityUtils.setSecurityManager(injector.getInstance(SecurityManager.class));
        }

        ClothoBuilder(Module... modulesArray){
            this(null, modulesArray);
        }
        
        private Injector injector;
        
        public <T> T get(Class<T> c){
            return injector.getInstance(c);
        }
        
        public <T> T get(Key<T> key){
            return injector.getInstance(key);
        }
    }
}
