 /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.HashMap;
import javax.script.AbstractScriptEngine;
import javax.script.Bindings;
import javax.script.Invocable;
import javax.script.ScriptContext;
import javax.script.ScriptEngineFactory;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import org.mozilla.javascript.*;
import static org.clothocad.core.execution.MetaEngine.API_NAME;

/**
 *
 * @author spaige
 */
class JavaScriptEngine extends AbstractScriptEngine implements HackEngine, Invocable {

    private final HashMap<Object, Object> indexedProps;
    private final ScriptableObject stdObjects;
    private final WrapFactory wrapFactory;

    public JavaScriptEngine() {
        //custom wrap factory
        wrapFactory = new WrapFactory() {
            @Override
            public Object wrap(Context cx, Scriptable scope, Object obj, Class<?> staticType) {
                final Object ret = super.wrap(cx, scope, obj, staticType);
                if (ret instanceof Scriptable) {
                    final Scriptable sret = (Scriptable) ret;
                    if (sret.getPrototype() == null) {
                        //set prototypes of Java Objects so that arbitrary
                        //properties on Objects can be set in the scripting environment
                        //they will not be persisted, though
                        sret.setPrototype(new NativeObject());
                    }
                }
                return ret;
            }
        };

        indexedProps = new HashMap<>();
        Context cx = Context.enter();
        cx.setWrapFactory(wrapFactory);
        stdObjects = cx.initStandardObjects();

        FunctionObject f = new JSLoader(this).getLoadFunction();
        stdObjects.put("load", stdObjects, f);
        cx.evaluateString(stdObjects, "var window = {};", "window hack", 1, null);
        cx.evaluateString(stdObjects, "var console = {}; console.log = function(){};", "console polyfill", 1, null);

        try {
            cx.evaluateReader(stdObjects, new FileReader("src/main/js/lib/lodash.js"), "lodash.js", 1, null);
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
            throw new RuntimeException("couldn't find lodash", ex);
        } catch (IOException ex) {
            throw new RuntimeException("couldn't open lodash.js", ex);
        } finally {
            Context.exit();
        }
    }

    @Override
    public Object eval(Reader reader, ScriptContext ctxt) throws ScriptException {
        Object ret;

        Context cx = Context.enter();
        cx.setWrapFactory(wrapFactory);
        try {
            Scriptable scope = getRuntimeScope(ctxt);
            String filename = "<Unknown source>";

            ret = cx.evaluateReader(scope, reader, filename, 1, null);
        } catch (RhinoException re) {
            re.printStackTrace();
            int line = (line = re.lineNumber()) == 0 ? -1 : line;
            String msg;
            if (re instanceof JavaScriptException) {
                msg = String.valueOf(((JavaScriptException) re).getValue());
            } else {
                msg = re.toString();
            }
            ScriptException se = new ScriptException(msg, re.sourceName(), line);
            se.initCause(re);
            throw se;
        } catch (IOException ee) {
            throw new ScriptException(ee);
        } finally {
            Context.exit();
        }

        return unwrapReturnValue(ret);
    }

    Scriptable getRuntimeScope(ScriptContext ctxt) {
        if (ctxt == null) {
            throw new NullPointerException("null script context");
        }

        // we create a scope for the given ScriptContext
        Scriptable newScope = new ExternalScriptable(ctxt, indexedProps);

        // Set the prototype of newScope to be 'topLevel' so that
        // JavaScript standard objects are visible from the scope.
        newScope.setPrototype(stdObjects);
        newScope.setParentScope(null);

        return newScope;
    }

    Object unwrapReturnValue(Object result) {
        if (result instanceof Wrapper) {
            result = ((Wrapper) result).unwrap();
        }

        return result instanceof Undefined ? null : result;
    }

    @Override
    public Object invokeFunction(String name, Object... args) throws ScriptException, NoSuchMethodException {
        return invoke(null, name, args);
    }

    @Override
    public Object invokeMethod(Object thiz, String name, Object... args) throws ScriptException, NoSuchMethodException {
        if (thiz == null) {
            throw new IllegalArgumentException("script object can not be null");
        }
        return invoke(thiz, name, args);
    }

    private Object invoke(Object thiz, String name, Object... args)
            throws ScriptException, NoSuchMethodException {
        Context cx = Context.enter();
        try {
            if (name == null) {
                throw new NullPointerException("method name is null");
            }

            if (thiz != null && !(thiz instanceof Scriptable)) {
                thiz = cx.toObject(thiz, stdObjects);
            }

            Scriptable engineScope = getRuntimeScope(context);
            Scriptable localScope = (thiz != null) ? (Scriptable) thiz
                    : engineScope;
            Object obj = ScriptableObject.getProperty(localScope, name);
            if (!(obj instanceof Function)) {
                throw new NoSuchMethodException("no such method: " + name);
            }

            Function func = (Function) obj;
            Scriptable scope = func.getParentScope();
            if (scope == null) {
                scope = engineScope;
            }
            Object result = func.call(cx, scope, localScope,
                    wrapArguments(args));
            return unwrapReturnValue(result);
        } catch (RhinoException re) {
            re.printStackTrace();
            int line = (line = re.lineNumber()) == 0 ? -1 : line;
            throw new ScriptException(re.toString(), re.sourceName(), line);
        } finally {
            cx.exit();
        }
    }

    Object[] wrapArguments(Object[] args) {
        if (args == null) {
            return Context.emptyArgs;
        }
        Object[] res = new Object[args.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = Context.javaToJS(args[i], stdObjects);
        }
        return res;
    }

    @Override
    public Bindings createBindings() {
        return new SimpleBindings();
    }

    @Override
    public ScriptEngineFactory getFactory() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Object eval(String script, ScriptContext context) throws ScriptException {
        return eval(new StringReader(script), context);
    }

    @Override
    public <T> T getInterface(Class<T> clasz) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public <T> T getInterface(Object thiz, Class<T> clasz) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void injectAPI(ScriptAPI api, ScriptContext context) {
        NativeJavaObject contextAPI = (NativeJavaObject) context.getAttribute(API_NAME);
        if (contextAPI != null) {
            // mutate preexisting api so we don't lose any custom properties
            ScriptAPI contextScriptAPI = (ScriptAPI) contextAPI.unwrap();
            contextScriptAPI.mind = api.mind;
            contextScriptAPI.api = api.api;
            //don't need to set persistor - persistor is singleton so will always be same instance
        } else {
            try{
                Context cx = Context.enter();
                cx.setWrapFactory(wrapFactory);
                context.setAttribute(API_NAME, wrapFactory.wrap(cx, getRuntimeScope(context), api, null), ScriptContext.ENGINE_SCOPE);
                cx.evaluateString(getRuntimeScope(context), "var require = function(string){return clotho.load(string)};", "alias clotho.load to require", 1, null);
            } finally {
                Context.exit();
            }
        }
    }
}
