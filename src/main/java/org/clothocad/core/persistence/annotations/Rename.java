/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.annotations;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 *
 * @author spaige
 */

@Retention(value = RetentionPolicy.RUNTIME)
public @interface Rename {
    public String value();
}
