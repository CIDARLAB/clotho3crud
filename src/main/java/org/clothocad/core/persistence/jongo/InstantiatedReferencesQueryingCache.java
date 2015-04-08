/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.BeanProperty;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.InjectableValues;
import com.fasterxml.jackson.databind.JavaType;
import com.fasterxml.jackson.databind.util.ClassUtil;
import com.mongodb.BasicDBObject;
import com.mongodb.DBCollection;
import com.mongodb.DBObject;
import com.thoughtworks.proxy.factory.CglibProxyFactory;
import com.thoughtworks.proxy.toys.hotswap.HotSwapping;
import com.thoughtworks.proxy.toys.hotswap.Swappable;
import java.util.HashSet;
import java.util.Set;
import javassist.Modifier;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.ReservedFieldNames;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.annotations.Reference;
import org.jongo.bson.Bson;
import org.jongo.marshall.MarshallingException;
import org.jongo.query.QueryFactory;

/**
 *
 * @author spaige
 */
@Slf4j
public class InstantiatedReferencesQueryingCache extends InjectableValues.Std{
    
    QueryFactory queryFactory;
    ExtendedUnmarshaller unmarshaller;
    DBCollection collection;
    
    public InstantiatedReferencesQueryingCache(QueryFactory queryFactory, ExtendedUnmarshaller unmarshaller, DBCollection collection){
        this.queryFactory = queryFactory;
        this.unmarshaller = unmarshaller;
        this.collection = collection;
    }
    
    public void addObjBase(ObjectId id, ObjBase value){
        addValue(id.toString(), value);
    }  
    
    @Override
    public Object findInjectableValue(Object valueId, DeserializationContext ctxt, BeanProperty forProperty, Object beanInstance) {
        //if a container, get content type
        
        JavaType type = forProperty.getType();
        if (type.isCollectionLikeType() | type.isArrayType()){
            type = type.getContentType();
        }
        //checking annotation catches properties with interface types
        if (ObjBase.class.isAssignableFrom(type.getRawClass()) || forProperty.getMember().getAnnotated().isAnnotationPresent(Reference.class)){                
            try {
                return super.findInjectableValue(valueId, ctxt, forProperty, beanInstance);
            } catch (IllegalArgumentException e){
                Object value = getValue(valueId, type);
                return value;
            }
            
        } else {
            return super.findInjectableValue(valueId, ctxt, forProperty, beanInstance);
        }
    }
    
    private ObjBase getValue(Object valueId, JavaType type){
        ObjectId id;
        try {
            id = new ObjectId(valueId);
        } catch (IllegalArgumentException e) {
            //TODO: collect these and similar warnings and send them to user output also
            log.warn("{} is an invalid ObjectId", valueId);
            return null;
        }
        DBObject result = query(id);
        if (result == null){
            log.warn("Could not find {} in database", id);
            return null;
        }        
        
        return makeValue(id, result, (Class<? extends ObjBase>) type.getRawClass());
    }

    
    public ObjBase makeValue(ObjectId id, DBObject rawData, Class<? extends ObjBase> c){

        if (Modifier.isAbstract(c.getModifiers()) || Modifier.isInterface(c.getModifiers())) {
            //extract declared schema from result, use that class
            Object schemaId = rawData.get(ReservedFieldNames.SCHEMA);
            try {
                c = (Class<ObjBase>) ClassUtil.findClass(schemaId.toString());
            } catch (ClassNotFoundException ex) {
                log.warn("Could not find fallback schema {}", schemaId);
                return null;
            }
        }
        
        ObjBase value = makeStub(id, c, rawData);
        // update stub with the rest of the data, 
        if (value != null && !(value instanceof Swappable)) 
            try {
                unmarshaller.unmarshall(Bson.createDocument(rawData), value, this);
            } catch (MarshallingException me){
                log.warn(me.getMessage(), me);
            }
        
        return value;
    }    
    
    public <T> T makeStub(ObjectId id, Class<T> type, DBObject rawData){
        T value;
        try {
            value = type.newInstance();
        } catch (IllegalAccessException | InstantiationException e) {
            try {
                value = makeWithJSONCreator(id, type, rawData);
            } catch (MarshallingException me) {
                log.warn(me.getMessage(), me);
                value = null;
            } catch (IllegalArgumentException iae){
                //make proxy to avoid circularity problems
                value = makeProxy(id, type);
            }
        }
        addValue(id.toString(), value);
        seen.remove(id);            
        return value;
    }

    @Override
    public Std addValue(String key, Object value) {
        //swap in real value if updating a proxy
        if (_values.containsKey(key) && _values.get(key) instanceof Swappable){
            Swappable proxy = (Swappable) _values.get(key);
            proxy.hotswap(value);
        }
        return super.addValue(key, value); //To change body of generated methods, choose Tools | Templates.
    }
    
    
    
    private DBObject query(ObjectId id){
        return collection.findOne(new BasicDBObject("_id", id.toString()));
    }
            
    private Set<ObjectId> seen = new HashSet<>();
    
    private <T> T makeWithJSONCreator(ObjectId id, Class<T> type, DBObject rawData){
        if (seen.contains(id)) 
            throw new IllegalArgumentException("Circularity when trying to construct " + id.toString());
        seen.add(id);
                //XXX: we end up unmarshalling the document twice in this case
        return unmarshaller.unmarshall(Bson.createDocument(rawData), type, this);
    }
    
    
    private <T> T makeProxy(ObjectId id, Class<T> c){
        T proxy = HotSwapping.proxy(c).with(null).build(new CglibProxyFactory());
        return proxy;
    }

}
