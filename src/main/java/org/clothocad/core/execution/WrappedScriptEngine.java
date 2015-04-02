/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import javax.script.Invocable;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptException;

/**
 *
 * @author spaige
 */
public class WrappedScriptEngine implements HackEngine {
    private final ScriptEngine engine;
    public WrappedScriptEngine(ScriptEngine engine){
        this.engine = engine;
    }

    @Override
    public Object eval(String script) throws ScriptException {
        return engine.eval(script);
    }
    
    @Override
    public Object invokeFunction(String name, Object... args) throws ScriptException, NoSuchMethodException {
        return ((Invocable) engine).invokeFunction(name, args);
    }

    @Override
    public Object invokeMethod(Object thiz, String name, Object... args) throws ScriptException, NoSuchMethodException {
        return ((Invocable) engine).invokeMethod(thiz, name, args);
    }

    @Override
    public Object eval(String script, ScriptContext context) throws ScriptException {
        return engine.eval(script, context);
    }
    
    @Override
    public void setContext(ScriptContext context){
        engine.setContext(context);
    }

    @Override
    public ScriptContext getContext() {
        return engine.getContext();
    }

    @Override
    public void injectAPI(ScriptAPI api, ScriptContext context) {
        if (engine instanceof HackEngine){
            ((HackEngine) engine).injectAPI(api, context);
        } else {
            context.setAttribute("clotho", api, ScriptContext.ENGINE_SCOPE);
        }
    }
}
