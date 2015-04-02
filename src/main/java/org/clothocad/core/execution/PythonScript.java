/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import lombok.Getter;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.IdUtils;

/**
 *
 * @author spaige
 */
public class PythonScript implements Script {
    //needs arg names, function name
    
    //how well does one engine for all functions scale?
    //do we need to manage the number of global functions? 
    
    //TODO: each script function needs to execute in its own scope
    //but
    //script functions should be cached somehow
    
    public PythonScript(){};
    
    public PythonScript(String source){
        this.source = source;
    }
    
    @Getter
    private String source;

    @Override
    public Set<ObjectId> findDependencies() {
       return new HashSet<>();
    }

    @Override
    public String generateImports(Collection<ObjectId> imports) {
        return "";
    }
    
    @Override 
    public String toString(){
        return source;
    }

    @Override
    public String modularizeFunction(String code) {
        return "";
    }

    
    @Override
    public String encapsulateModule(String code, String setupcode) {
        return "";
    }
}
