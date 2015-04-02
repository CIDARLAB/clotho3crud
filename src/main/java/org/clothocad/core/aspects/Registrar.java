/*
 * 
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.aspects;

import java.util.Map;
import org.clothocad.core.persistence.Persistor;
/**
 * The Registrar (which maybe should be renamed) is the Clotho-to-Clotho API
 * 
 * @author John Christopher Anderson
 */


public class Registrar {
    /**
     * Clotho-to-Clotho request
     * for confirmation that a given user has been given a particular badge
     * 
     * @param personId
     * @param badgeId
     * @return 
     */
    
    Persistor persistor;
    
    public boolean hasBadge(String personId, String badgeId) {
        try {
        	return true;
        	/** TODO!!
            Badge badge = persistor.get(Badge.class, new ObjectId(badgeId));
            return badge.hasBadge(personId);
            **/
        } catch(Exception err) {
            return false;
        }
    }
    
    
    //HMMM, THIS IS BEGINNING TO FEEL REDUDANT WITH THE SERVERSIDE API, THINGS SHOULD BE RELAYED ONTO AN INSTANCE OF THAT I THINK
    
    //MAYBE ANOTHER INSTANCE OF CLOTHO IS SIMPLY ISSUED A MIND PLUS SOME DIRECT JAVA-LEVEL ADDITIONAL METHODS FOR AUTHENTICATION?
    
    /**
     * Clotho=to-Clotho request
     * for a Sharable object that exists on this domain
     * 
     * @param personId
     * @param sharableId
     * @return 
     */
    public Map<String, Object> get(String personId, String sharableId) {
        try {
        	/*** TODO!
            Sharable out = persistor.get(Sharable.class, new ObjectId(sharableId));
//            if(out.canShare(personId)) {
            if(true) {
                return out.toJSON();
            }
            ***/
            //return null;
        } catch(Exception err) {
        }
        return null;
    }
}
