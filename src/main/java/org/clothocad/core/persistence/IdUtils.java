/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import java.util.Map;
import javax.inject.Inject;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.schema.BuiltInSchema;
import org.clothocad.core.schema.Schema;

/**
 * lame hack to handle the need to look up various ids in the data layer
 * persistor injected by clotho module
 * @author spaige
 */
public class IdUtils {
    @Inject
    static Persistor persistor;
    
    public static ObjBase get(ObjectId id){
        return persistor.get(ObjBase.class, id);
    }
    
    public static Class getClass(ObjectId schemaId) throws ClassNotFoundException{
        Map<String,Object> schemaData = persistor.getAsJSON(schemaId);
        if (schemaData.get("schema").toString().equals("BuiltInSchema")){
            BuiltInSchema schema = persistor.get(BuiltInSchema.class, schemaId);
            return schema.getClass();
        }
        else {
            return Class.forName(schemaId.toString(), true, Schema.cl);
        }
    }
}
