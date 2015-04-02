/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.fasterxml.jackson.annotation.JsonView;
import java.util.Map;
import java.util.Set;
import javax.inject.Inject;
import lombok.Data;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.persistence.annotations.Add;
import org.clothocad.core.persistence.annotations.Adds;
import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.jackson.JSONViews;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.util.JSON;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 */

@EqualsAndHashCode(exclude={"fields", "methods"}, callSuper = true)
@Data
@NoArgsConstructor
@Slf4j
@Adds({@Add(name="language", provider="getLanguage")})
public abstract class Schema extends SharableObjBase {
    
    public Schema(String name, String description, Person author){
        super(name, author, description);
    }
    
    protected static final String BASE_PACKAGE_BINARY = "org.clothocad.loadedschemas.";
   
    @Inject
    public static  DBClassLoader cl = null;
    
    @JsonView(JSONViews.Internal.class)
    protected byte[] classData;
    protected Map<String, ObjectId> dependencies;
    protected String source;
    
    //These are settable only in ClothoSchema - they are derived from source in other languages
    
    protected Set<ClothoField> fields;
    protected Set<Function> methods;
    
    @Reference
    protected Schema superClass;

    public abstract Language getLanguage();

    //TODO: handle files that result in multiple source files;
    //TODO: supply uuid's of referenced classes to resolve name conflicts
    public abstract void setSource(String source);
    //***Proxying methods to method handles???
    //can get bytecode from functions? 
   
    public String getBinaryName(){
        return this.getId().toString();
    }
    
    public String getInternalName(){
        return getBinaryName().replace('.', '/');
    }
    
    //the classloader can only find saved schemas, so if this throws an exception, try saving the schema
    public Class<? extends ObjBase> getEnclosedClass(ClassLoader cl) throws ClassNotFoundException{
        return (Class<? extends ObjBase>) cl.loadClass(getBinaryName());
    }

    public boolean childOf(Schema schema) {
        if (schema == null || superClass == null) return false;
        if (schema.getId().equals(superClass.getId())) return true;
        return (childOf(schema.superClass));
    }
    
    public ObjBase instantiate(Map<String,Object> data) throws ClassNotFoundException{
        Class<? extends ObjBase> c = getEnclosedClass(cl);
        return JSON.mapper.convertValue(data, c);
    }
}
