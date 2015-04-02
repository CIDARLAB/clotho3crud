/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import com.fasterxml.jackson.annotation.JsonIgnore;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.reflect.Member;
import javax.script.ScriptException;
import org.mozilla.javascript.Context;
import org.mozilla.javascript.FunctionObject;
import org.mozilla.javascript.Scriptable;
import org.mozilla.javascript.ScriptableObject;

/**
 *
 * @author spaige
 */
public class JSLoader extends ScriptableObject {

    public JSLoader(JavaScriptEngine engine) {
        this.engine = engine;
    }
    
    JavaScriptEngine engine;
    
    @Override
    public String getClassName() {
        return getClass().getName();
    }

    //XXX:
    //There are a lot of problems with this approach
    // pollutes global namespace
    // doesn't protect against multiple imports
    // we should move to something like rhino-require
    public void load(String filename) throws ScriptException{
        try {
            engine.eval(new FileReader("clotho3-web/lib/"+filename));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Error loading javascript file: " + filename, e);
        }
    }
    
    @JsonIgnore
    public FunctionObject getLoadFunction(){
        try {
            return new InstanceBoundFunctionObject("load", JSLoader.class.getDeclaredMethod("load", String.class), this);
        } catch (NoSuchMethodException ex) {
            ex.printStackTrace();
        } catch (SecurityException ex) {
            ex.printStackTrace();
        }
        return null;
    }
    
    
    //http://stackoverflow.com/questions/3441947/how-do-i-call-a-method-of-a-java-instance-from-javascript
    
    private static class InstanceBoundFunctionObject extends FunctionObject {

        private InstanceBoundFunctionObject(String name, Member methodOrConstructor, Scriptable parentScope) {
            super(name, methodOrConstructor, parentScope);
        }

        @Override
        public Object call(Context cx, Scriptable scope, Scriptable thisObj, Object[] args) {
            return super.call(cx, scope, getParentScope(), args);
        }
    }
}
