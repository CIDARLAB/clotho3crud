/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import com.fasterxml.jackson.annotation.JsonAutoDetect;
import com.fasterxml.jackson.annotation.JsonIgnore;
import java.util.HashSet;
import java.util.Set;
import lombok.Getter;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.util.Language;
import static org.clothocad.core.datums.util.Language.JAVASCRIPT;
import org.clothocad.core.execution.JavaScriptScript;
import org.clothocad.core.execution.PythonScript;
import org.clothocad.core.execution.Script;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

/**
 *
 * @author spaige
 */
@JsonAutoDetect(getterVisibility = JsonAutoDetect.Visibility.PUBLIC_ONLY, 
                setterVisibility = JsonAutoDetect.Visibility.PUBLIC_ONLY)
@Slf4j
public class Module extends ObjBase {
    
    protected Script code;
    @ReferenceCollection
    protected Module[] dependencies;
    protected String description;
    @Getter
    //it's currently unclear to me whether modules or scripts or both 
    //should have the language parameter
    
    //right now scripts are just serialized to the string representation of their source code, 
    // which is nice and simple
    // but that means we can't have modules with different languages inside them
    protected Language language;

    public Module() {
    }
    public Module(String name, String description, Language language, 
            String code, Module[] dependencies){
        setName(name);
        this.description = description;
        this.language = language;
        this.dependencies = dependencies;
        setCode(code);
    }

    public void setLanguage(Language language){
        this.language = language;
        //check for cached code
        if (cachedCode != null){
            setCode(cachedCode, language);
            cachedCode = null;
        }
    }
    
    private transient String cachedCode;
    
    public void setCode(String code){
        if (language != null) setCode(code, language);
        else {
            //cache to set later when language is set
            cachedCode = code;
        }
    }
    
    protected void setCode(String code, Language language){
        this.code = generateScript(code,language);
    }
    
    public static Script generateScript(String code, Language language){
        switch (language){
            case JAVASCRIPT:
                return new JavaScriptScript(code);
            case PYTHON:
                return new PythonScript(code);
            default:
                throw new UnsupportedOperationException("unsupported language");
        }          
    }

    //TODO: test w/ @PostLoad also
    protected void syncDependencies(){
        //figure out dependencies declared in code 
        Set<ObjectId> declaredDependencies = code.findDependencies();
        Set<ObjectId> listedDependencies = getDependencySet(dependencies);
      
        Set<ObjectId> listedButNotDeclared = new HashSet<>(listedDependencies);
        listedButNotDeclared.removeAll(declaredDependencies);
        //Set<ObjectId> declaredButNotListed;
        
        //declare all listed dependencies in code
        //code.addImports(listedButNotDeclared);
        //TODO: add declared but not listed dependencies to dependency list
    }
    
    private static Set<ObjectId> getDependencySet(Module[] dependencies) {
        Set<ObjectId> output = new HashSet<>();
        
        if (dependencies != null) for (Module obj : dependencies){
            output.add(obj.getId());
        }
        return output;
    }
    
    public String getCode() {
        if (language == null) return cachedCode;
        return code == null ? null : code.toString();
    }
    
    @JsonIgnore
    public String getCodeToLoad() {
        return code.encapsulateModule(code.getSource(), getSetup());
    }
    
    @JsonIgnore
    public String getSetup(){
        syncDependencies();
        return this.code.generateImports(getDependencySet(dependencies));
    }

    @JsonIgnore
    public Function getFunction(String name) {
        Function function = new Function(name, null, null, null, language);
        function.dependencies = new Module[]{this};
        String formatString = "function () { return %s.%s.apply(this,arguments);};";
        function.setCode(String.format(formatString, this.getName(), name));
        return function;
    }
}
