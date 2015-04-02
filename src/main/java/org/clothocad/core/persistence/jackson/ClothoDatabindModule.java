/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.core.Version;
import com.fasterxml.jackson.databind.BeanDescription;
import com.fasterxml.jackson.databind.DeserializationConfig;
import com.fasterxml.jackson.databind.Module;
import com.fasterxml.jackson.databind.deser.BeanDeserializerBuilder;
import com.fasterxml.jackson.databind.deser.BeanDeserializerModifier;

/**
 *
 * @author spaige
 */
public class ClothoDatabindModule extends Module{

    private static final Version version = new Version(0,0,0, "", "org.clothocad", "clotho-databind");
    
    @Override
    public String getModuleName() {
        return "Clotho Databind";
    }

    @Override
    public Version version() {
        return version;
    }

    @Override
    public void setupModule(SetupContext context) {
        context.addBeanDeserializerModifier(new BeanDeserializerModifier(){
            @Override
            public BeanDeserializerBuilder updateBuilder(DeserializationConfig config, BeanDescription beanDesc, BeanDeserializerBuilder builder) {
                //TODO: use clotho Value Instantiator (gets objects from pre-processing)
                builder.setValueInstantiator(null);
                return builder;
            }
        });
    }
    
}
