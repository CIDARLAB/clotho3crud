/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import javax.script.ScriptContext;
import javax.script.ScriptException;

/**
 *
 * @author spaige
 */
public interface HackEngine {
    public Object eval(String script) throws ScriptException;
    public Object eval(String script, ScriptContext context) throws ScriptException;
    
    public Object invokeFunction(String name, Object... args)
        throws ScriptException, NoSuchMethodException;
    
    public Object invokeMethod(Object thiz, String name, Object... args)
        throws ScriptException, NoSuchMethodException; 
    
    public void setContext(ScriptContext context);
    public ScriptContext getContext();
    
    public void injectAPI(ScriptAPI api, ScriptContext context);
}
