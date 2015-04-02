/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import org.clothocad.model.Person;

/**
 *
 * @author jcanderson
 */
public interface Sharable  {
    //Metadata for all Sharables
    public ObjectId getId();
    public Person getAuthor();
    public String getIcon();
    public String getName();
    public String getDescription();
    
}