/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

//import com.sun.script.javascript.RhinoScriptEngine;
import javax.script.ScriptContext;
//import javax.script.ScriptEngine;
//import lombok.Delegate;
//import org.mozilla.javascript.Context;
import org.mozilla.javascript.Scriptable;

/**
 *
 * @author spaige
 * 
 * Temp fix: Large sections commented out b/c build failing, compiler cannot find
 * com.sun.script.javascript
 */

public class ClothoRhinoScriptEngine { //implements ScriptEngine {
    
    //@Delegate
    //private RhinoScriptEngine engine;
    
    Scriptable getRuntimeScope(ScriptContext ctxt) {
        if (ctxt == null) {
            throw new NullPointerException("null script context");
        }
        Scriptable newScope = null;
        // we create a scope for the given ScriptContext
        //TODO: introspection shenanigans Scriptable newScope = new ExternalScriptable(ctxt, indexedProps);

        // Set the prototype of newScope to be 'topLevel' so that
        // JavaScript standard objects are visible from the scope.
        //TODO: introspection shenanigans newScope.setPrototype(topLevel);

        // define "context" variable in the new scope
        newScope.put("context", newScope, ctxt);

        //TODO?: define print, println
        return newScope;
    }
   
}
