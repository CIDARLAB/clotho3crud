/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import com.google.common.base.Joiner;
import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import org.clothocad.core.datums.util.Language;


public class CompiledEngineScript extends ScriptEngineScript {

    public CompiledEngineScript() {
    }

    public CompiledEngineScript(String name,  String[] argNames, String source, Language language) {
        super(name, source, language);
        //this.returnValue = returnValue;
        //this.argNames = argNames;
        this.argNames = argNames;
    }
    
    private String returnValue;
    private String[] argNames;
    private CompiledScript compiled;
    private Bindings env;
    
    @Override
    protected void load() throws ScriptException{
        env = new SimpleBindings();
        engine.eval(source, env);
        
        
        compiled = ((Compilable) engine).compile(name + "(" + Joiner.on(",").join(argNames) + ")");
        loaded = true;
        //TODO: dependencies inserted into Bindings.
    }
    
    @Override
    public Object run(Object... args) throws ScriptException{
        //todo: Throw exception if args.lengh != argNames.length
        if (!loaded) load ();
        for (int i =0;i< args.length;i++){
            env.put(argNames[i], args[i]);
        }
        Object result = compiled.eval(env);
        return result;
    }
}
