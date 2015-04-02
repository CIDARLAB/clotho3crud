/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author spaige
 */
public class ReflectionUtils {
    public static Method findMethodNamed(String name, int argsLength, Class type){
        for (Method method : type.getMethods()){
            if (method.getName().equals(name) && method.getParameterTypes().length == argsLength){
                return method;
            }
        }
        return null;
    }
    
    public static Map<String, Object> getFieldsAndValues(Object o){
        Map<String,Object> out = new HashMap<>();
        Field[] fields = o.getClass().getFields();
        for (Field field : fields){
            try {
                out.put(field.getName(), field.get(o));
            } catch (IllegalArgumentException ex) {
                ex.printStackTrace();
            } catch (IllegalAccessException ex) {
                ex.printStackTrace();
            }
        }
        return out;
    }
}
