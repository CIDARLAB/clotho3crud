/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineFactory;
import javax.script.ScriptEngineManager;    
import org.clothocad.core.datums.util.Language;

/**
 *
 * @author spaige
 */
public class ClothoScriptEngineManager {
    private static ScriptEngineManager em = new ScriptEngineManager();
    
    private static Map<Language, ScriptEngine> enginesByLanguage = new HashMap<>();
    
    static {
        Properties props = System.getProperties();
        //props.setProperty("rhino.opt.level", "9");
    }

    private static ScriptEngineFactory getFactoryByLanguage(Language language){
        for (ScriptEngineFactory factory : em.getEngineFactories()){
            for (String extension : factory.getExtensions()){
                if (language.extensions().contains(extension)){
                    return factory;
                }
            }
        }   
        return null;
    }
        
    public static ScriptEngine getEngineByLanguage(Language language){
        if (language == Language.JAVASCRIPT) return new JavaScriptEngine();

        ScriptEngineFactory factory = getFactoryByLanguage(language);
        boolean threadSafe = (factory.getParameter("THREADING") != null);
        if (threadSafe && enginesByLanguage.containsKey(language)){
            return enginesByLanguage.get(language);
        }
        ScriptEngine engine = factory.getScriptEngine();
        
        if (threadSafe) enginesByLanguage.put(language, engine);
        return engine;
    }
}
