/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 */
@NoArgsConstructor
public abstract class SharableObjBase extends ObjBase implements Sharable {
    
    public SharableObjBase(String name, Person author){
        setName(name);
        this.author = author;
    }

    public SharableObjBase(String name, Person author, String description){
        this(name, author);
        this.description = description;
    }
    
    @Getter @Setter
    private Person author;

    @Getter @Setter
    private String description, icon;
    
}
