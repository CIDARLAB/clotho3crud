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

package org.clothocad.core.datums;

import java.util.HashSet;
import java.util.Set;

/**
 * This the serverside implementation of the Client
 * Both the Client and an associated User must be logged in to connect
 * A client is at any time associated with one user

 * @author John Christopher Anderson
 */


public class Client 
		extends ObjBase {
    
    private String macAddress; //The actual UUID of the client
    private String userId;
    private boolean isLoggedIn = false; //The maintainence of login state for this Client
    private Set<User> users = new HashSet<User>(); //Stores the links to the users that have used this Client
            
    void lockout() {
        isLoggedIn = false;
        
        //Persistor.get().persist(this);
        // new:
        //this.save();
    }
    
    void login() {
        isLoggedIn = true;
        //Persistor.get().persist(this);
        // new:
        //this.save();
    }
}
