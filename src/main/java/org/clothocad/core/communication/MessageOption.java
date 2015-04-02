/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

/**
 *
 * @author spaige
 */
public enum MessageOption {
    
    detail, //ID_ONLY (returns only the id), NORMAL (returns object description)
    mute, //true, 
    maxResults, //maximum results number
    filter; //list of field names to filter on
    
    public final static String ID_ONLY = "ID_ONLY";
    public final static String NORMAL = "NORMAL";
}
