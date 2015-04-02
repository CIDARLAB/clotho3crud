/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import javax.script.ScriptException;
import org.clothocad.core.datums.Function;

/**
 *
 * @author spaige
 */
public class JavaScriptAPI extends AbstractScriptAPI {
    
    @Override
    public void injectFunction(Function f){
        //TODO: name conflict detection?
        
        //figure out how to put a Rhino function in the scripter that proxies execute instead of injecting the object and a function to call it
        engine.put(f.getName(), f);
        try {
            //if you do this inside a method call, in which scope does the function end up?
            engine.eval(f.getName() + " = " + f.getName()+".execute");
        } catch (ScriptException ex) {
            throw new RuntimeException(ex);
        }
    }
    
}
