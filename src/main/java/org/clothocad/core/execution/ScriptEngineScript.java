/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import java.util.Collection;
import java.util.Set;
import javax.script.Invocable;
import javax.script.ScriptEngine;
import javax.script.ScriptException;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.util.Language;

/**
 *
 * @author spaige
 */

public class ScriptEngineScript implements Script {
        //needs arg names, function name
    
    //how well does one engine for all functions scale?
    //do we need to manage the number of global functions? 
    
    //TODO: each script function needs to execute in its own scope
    //but
    //script functions should be cached somehow
    protected transient ScriptEngine engine;
       
    protected String source;
    protected String name;
    private Language language;
    
    public ScriptEngineScript(){};
    
    public ScriptEngineScript(String name, String source, Language language){
        this.name = name;
        this.source = source;
        this.language = language;
        
        getScriptEngine();
    }
    
    protected transient boolean loaded = false;
    
    protected void load() throws ScriptException{
        Object result = null;
        try {
            result = engine.get(name);
        } catch(IllegalArgumentException e){
        }
        if (result == null){
            engine.eval(source);
        }
        
        loaded = true;
    }
    
    //XXX: replace @PostLoad w/ non-morphia annotation
    public void getScriptEngine(){
        engine = ClothoScriptEngineManager.getEngineByLanguage(language);
    }
    
    
    
    public Object run(Object... args) throws ScriptException {
        Invocable invocable = (Invocable) engine;
        if (!loaded) load();
        try {
            return invocable.invokeFunction(name, args);
        } catch (NoSuchMethodException ex) {
            throw new RuntimeException(ex);
        }
    }

    @Override
    public Set<ObjectId> findDependencies() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String getSource() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String generateImports(Collection<ObjectId> listedButNotDeclared) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String modularizeFunction(String code) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String encapsulateModule(String code, String setupCode) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
