/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import java.util.concurrent.Callable;

/**
 *
 * @author spaige
 */
public class InitializePersistor implements Callable {
    private Persistor persistor;
    
    public InitializePersistor(Persistor persistor){
        this.persistor = persistor;
    }
    @Override
    public Object call() throws Exception {
        persistor.initializeBuiltInSchemas();
        return null;
    }
     
}
