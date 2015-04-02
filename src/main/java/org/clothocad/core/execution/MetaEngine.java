/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.annotation.JsonSerialize;
import com.fasterxml.jackson.databind.util.StdConverter;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.script.Bindings;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import javax.script.SimpleScriptContext;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Module;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.ObjectId;
import org.mozilla.javascript.RhinoException;

/**
 *
 * @author spaige
 */
public class MetaEngine {
    public static final String API_NAME = "clotho";

    private transient Map<Language, ScriptContext> contexts = new HashMap<>();
    private transient Map<Language, ScriptEngine> engines = new HashMap<>();
    
    @JsonSerialize(converter = BindingsMapConverter.class)
    @JsonDeserialize(contentAs = SimpleBindings.class)
    private Map<Language, Bindings> bindings = new HashMap<>();

    public static class BindingsMapConverter extends StdConverter<Map<Language,Bindings>,Map<Language,Bindings>>{

        @Override
        public Map<Language,Bindings> convert(Map<Language,Bindings> value) {
            for (Language key : value.keySet()){
                Bindings bindings = value.get(key);
                //api does not need to be serialized, and does not serialize cleanly
                bindings.remove("clotho");
                bindings.remove("require");
            }
            return value;
        }
        
    }
    
    public Object eval(String script, Language language, ScriptAPI api) throws ScriptException {
        ScriptContext context = getContext(language);
        HackEngine engine = getEngine(language);

        engine.injectAPI(api, context);
        try{
            return engine.eval(script, context);
        } //XXX: ew, referring to a concrete engine implementation!
        catch (RhinoException e) {
            //is this terrible?
            throw new ScriptException(e);
        }
    }

    private ScriptContext getContext(Language language) {
        if (!contexts.containsKey(language)) {
            ScriptContext context = new SimpleScriptContext();
            if (bindings.containsKey(language)) {
                context.setBindings(bindings.get(language), ScriptContext.ENGINE_SCOPE);
            } else {
                bindings.put(language, context.getBindings(ScriptContext.ENGINE_SCOPE));
            }
            contexts.put(language, context);
        }
        return contexts.get(language);
    }

    private HackEngine getEngine(Language language) {

        if (!engines.containsKey(language)) {
            engines.put(language, ClothoScriptEngineManager.getEngineByLanguage(language));
            ScriptEngine engine = engines.get(language);
            if (engine == null) {
                throw new EngineNotFoundException();
            }
            engine.setContext(getContext(language));
        }
        return new WrappedScriptEngine(engines.get(language));
    }

    public Object get(String token) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    public Object invoke(String functionCode, String name, List args, ScriptAPI api) throws ScriptException {
        return invoke(functionCode, "", name, args, api);
    }

    public static String idToName(ObjectId id) {
        return "f" + id.toString().replaceAll("\\W", "_");
    }

    public void loadAsGlobal(Module module, ScriptAPI api) throws ScriptException {
        String name = idToName(module.getId());
        //TODO: check last saved date
        //if (getEngine(module.getLanguage()).getContext().getBindings(ScriptContext.ENGINE_SCOPE).containsKey(name)){
        //    return;
        //}
        HackEngine engine = getEngine(module.getLanguage());

        try {
            engine.eval("var " + name + " = " + module.getCodeToLoad(), engine.getContext());
        } catch (ScriptException e){
            throw new ScriptException(e.getMessage(),module.getId().toString(),e.getLineNumber(), e.getColumnNumber());
        }
    }

    public Object invoke(Function function, List args, ScriptAPI api) throws ScriptException, NoSuchMethodException {
        HackEngine engine = getEngine(function.getLanguage());
            engine.injectAPI(api, engine.getContext());
            loadAsGlobal(function, api);
            return engine.invokeFunction(idToName(function.getId()), args.toArray());
    }

    public Object invoke(Module module, String methodName, List args, ScriptAPI api) throws ScriptException, NoSuchMethodException {
        HackEngine engine = getEngine(module.getLanguage());
       engine.injectAPI(api, engine.getContext());
            loadAsGlobal(module, api);

            Object thiz = engine.getContext().getBindings(ScriptContext.ENGINE_SCOPE).get(idToName(module.getId()));
            return engine.invokeMethod(thiz, methodName, args.toArray());
    }

    //XXX: de-js-ify this!
    public Object invoke(String functionCode, String setupCode, String name, List args, ScriptAPI api) throws ScriptException {
       
            HackEngine engine = getEngine(Language.JAVASCRIPT);
            ScriptContext context = engine.getContext();
            engine.injectAPI(api,context);

            //eval dependencies
            //XXX: this is terrible terrible, change to module pattern that returns a function
         try {
            engine.eval(setupCode);

            engine.eval("var " + name+ " = " + functionCode);
                
            
            return engine.invokeFunction(name, args.toArray());

        } catch (NoSuchMethodException ex) {
            throw new ScriptException(ex);
        } 
    }

    protected String generateDependencies(Language language, Collection<ObjectId> dependencies) {
        switch (language) {
            case JAVASCRIPT:
                return new JavaScriptScript().generateImports(dependencies);
            default:
                throw new UnsupportedOperationException();
        }
    }
}
