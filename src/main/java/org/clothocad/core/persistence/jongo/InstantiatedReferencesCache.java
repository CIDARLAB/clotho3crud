/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.BeanProperty;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.InjectableValues;
import com.fasterxml.jackson.databind.JavaType;
import com.thoughtworks.proxy.factory.CglibProxyFactory;
import com.thoughtworks.proxy.toys.hotswap.HotSwapping;
import com.thoughtworks.proxy.toys.hotswap.Swappable;
import java.util.ArrayList;
import java.util.List;
import javassist.Modifier;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.util.BasicObjBase;

/**
 *
 * @author spaige
 */
public class InstantiatedReferencesCache extends InjectableValues.Std{
    
    List<ObjBase> undone = new ArrayList<>();
    
    public void addObjBase(String valueId, ObjBase value){
        addValue(valueId, value);
        undone.add(value);
    }  
    
    public void addProxy(String valueId, ObjBase value){
        addValue(valueId, value);
        undone.add(0, value);
    }
    
    public void resolveProxy(String valueId, ObjBase value){
        try {
            Swappable proxy = Swappable.class.cast(_values.get(valueId));
            proxy.hotswap(value);
        } catch (ClassCastException e) {
            throw new RuntimeException("Internal Error: tried to resolve a non-proxy object", e);
        }
    }

    @Override
    public Object findInjectableValue(Object valueId, DeserializationContext ctxt, BeanProperty forProperty, Object beanInstance) {
        //if a container, get content type
        
        JavaType type = forProperty.getType();
        if (type.isCollectionLikeType() | type.isArrayType()){
            type = type.getContentType();
        }
        if (ObjBase.class.isAssignableFrom(type.getRawClass())){
            Class<ObjBase> c = (Class<ObjBase>) type.getRawClass();
            try {
                //already made an object for this id - return it
                return super.findInjectableValue(valueId, ctxt, forProperty, beanInstance);
            } catch (IllegalArgumentException e){
                // else, create it, cache it, return it
                //XXX: just assuming there's a lombok no-arg constructor for now
                if (Modifier.isAbstract(c.getModifiers()) || Modifier.isInterface(c.getModifiers())){
                    //todo: also if no-args
                    
                    //make a proxy and return that
                    ObjBase newValue = HotSwapping.proxy(c).with(new BasicObjBase()).build(new CglibProxyFactory());
                    this.addProxy(valueId.toString(), newValue);
                    newValue.setId(new ObjectId(valueId.toString()));
                    return newValue;
                }
                
                try{
                    ObjBase newValue = c.newInstance();
                    newValue.setId(new ObjectId(valueId));
                    this.addObjBase(valueId.toString(), newValue);
                    return newValue;
                } catch (InstantiationException|IllegalAccessException ie){
                    throw new RuntimeException("Could not invoke no-args constructor on " + c.getCanonicalName(), ie);
                }
            }
            
        } else {
            return super.findInjectableValue(valueId, ctxt, forProperty, beanInstance);
        }
    }
    
    public Boolean done(){
        return undone.isEmpty();
    }
    
    public ObjBase getNextUndone(){
        return undone.remove(0);
    }

    boolean has(String id) {
        return _values.containsKey(id);
    }
    
    public Object get(String id){
        return _values.get(id);
    }
}
