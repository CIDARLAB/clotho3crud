/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.schema;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 *
 * @author prashantvaidyanathan
 */
public class RunScripts {
    
    public Object runningfunction(Object obj,Method func, String message) throws Exception
    {
        Object[] parameters = new Object[1];
        parameters[0] = message;
        return func.invoke(obj, parameters);
    }
    
    public Object runFunction(Object classobj, Method method, Object[] arguements) throws Exception
    {
        /*Object[] parameters;
        parameters = new Object[arguements.length];
        
        for(int i=0;i<parameters.length;i++)
        {
            parameters[i] = arguements[i];
        }*/
        return method.invoke(classobj, arguements);
    }
}

