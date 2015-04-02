/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.util.Language;

/**
 * For schemas that don't actually exist, but have been referred to by name by 
 * saved data
 *
 * @author spaige
 */


public class InferredSchema extends Schema{

    public InferredSchema(String name){
        this.setName(name);
    }
    

    @Override
    public Language getLanguage() {
        return Language.JSONSCHEMA; //not quite accurate
    }

    @Override
    public void setSource(String source) {
        //TODO: replaces self with real schema?
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getBinaryName() {
        return null;
    }
    
    @Override
    public String getInternalName(){
        return null;
    }

    @Override
    public Class<? extends ObjBase> getEnclosedClass(ClassLoader cl) throws ClassNotFoundException {
        return null;
    }
    
    
}
