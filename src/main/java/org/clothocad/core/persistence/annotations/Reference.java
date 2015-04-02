/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.annotations;

import com.fasterxml.jackson.annotation.JacksonAnnotationsInside;
import com.fasterxml.jackson.annotation.JsonTypeInfo;
import com.fasterxml.jackson.annotation.JsonTypeInfo.Id;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.annotation.JsonSerialize;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import org.clothocad.core.persistence.jackson.ReferenceDeserializer;
import org.clothocad.core.persistence.jackson.ReferenceSerializer;

/**
 * Annotation flag marking a property that should be serialized/deserialized 
 * as reference to another entry, or a property that is a container of such values
 * 
 * Contains jackson annotations for directing ObjectMapper to appropriate 
 * serializer/deserializer
 * 
 * @author spaige
 */
@Retention(RetentionPolicy.RUNTIME)
@JacksonAnnotationsInside
@JsonSerialize(using=ReferenceSerializer.class)
@JsonDeserialize(using=ReferenceDeserializer.class)
@JsonTypeInfo(use = Id.NONE)
public @interface Reference {
    
}
