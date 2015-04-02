/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.annotations;

import java.lang.annotation.ElementType;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 *
 * @author spaige
 */
@Target(value = {ElementType.TYPE})
@Retention(value = RetentionPolicy.RUNTIME)
@Inherited
public @interface Add {
    public String name();
  
    //XXX: not supported yet
    public String value() default "";
    
    public String provider() default "";
    
    public boolean isReference() default false;
    
    public Class concreteClass() default Object.class;
    
    //XXX: doesn't do anything until change tracking is implemented
    public String[] dependsOn() default {};
}
