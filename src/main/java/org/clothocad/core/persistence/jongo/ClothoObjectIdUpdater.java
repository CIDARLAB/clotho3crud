/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import java.lang.reflect.Field;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.jongo.ReflectiveObjectIdUpdater;
import org.jongo.marshall.jackson.JacksonIdFieldSelector;

/**
 *
 * @author spaige
 */
public class ClothoObjectIdUpdater extends ReflectiveObjectIdUpdater{

    public ClothoObjectIdUpdater(IdFieldSelector idFieldSelector) {
        super(idFieldSelector);
    }

    @Override
    public Object getId(Object pojo) {
        if (pojo instanceof ObjBase){
            ObjBase obj = (ObjBase) pojo;
            ObjectId id = obj.getId();
            return obj.getId().toString();
        }
        return super.getId(pojo); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setObjectId(Object newPojo, org.bson.types.ObjectId id) {
        if (newPojo instanceof ObjBase){
            ObjBase obj = (ObjBase) newPojo;
            obj.setId(new ObjectId(id.toString()));
            return;
        }
        super.setObjectId(newPojo, id); //To change body of generated methods, choose Tools | Templates.
    }
    
    public static class ClothoIdFieldSelector extends JacksonIdFieldSelector{

        @Override
        public boolean isId(Field f) {
            //it would be _better_ if we checked for mixins, but that's not feasible right now
            if (ObjBase.class.isAssignableFrom(f.getDeclaringClass())){
                return (f.getName().equals(ID));
            }
            return super.isId(f); //To change body of generated methods, choose Tools | Templates.
        }
        
    }
    
}
