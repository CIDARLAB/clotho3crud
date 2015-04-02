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
import java.util.UUID;
import org.clothocad.core.persistence.Persistor;

/**
 * User's are used to authenticate individuals (by the Authenticator)
 * It's the object that represents an individual and whether or not
 * they are globally logged off or not.
 * 
 * Note that User explicitly does not implement Sharable.
 * 
 * The User is also entirely inaccessible from within the container.
 * 
 * When you log in to Clotho you are loggin in as a person, which
 * is relayed onto a user.  This is a personal security measure to allow people
 * to make distinct separations between modes of use.  For example, if an
 * individual works for two employers, and their IP needs to stay separated,
 * then that user would create two separate Persons to log in with.
 * 
 * When logged in, all aspects of your world are related to your person, not
 * your user.  The core does know about your Persons via their uuid, but the container
 * must be indexed from within the container, so indexing will occur over Persons,
 * not Users.
 * 
 * I think the deal is that Persons are normal dynamic objects referred to by the core
 * It's the one piece of the container hard coded by the core.  Other than that it is 
 * a normal Schema with normal ObjBases that follow normal Sharable rules.  So, Person
 * objects pass freely across domains (if you allow that).  However, Users own objects
 * while Persons have permissions for objects and are the owners of Objects (since people
 * will need to see who made stuff from within the container).
 * For security reasons, a User should never leak out of its domain.  You create a User on
 * a per-domain basis, and on each domain your Person maps to a User.  So, login state is 
 * referential to a particular domain.
 *
 * @author John Christopher Anderson
 */


public class User 
	extends ObjBase {
	
    private String authKey;
    
    public User() {
        //YEAH, I DON'T KNOW WHAT THIS AUTHENTICATION KEY IS ABOUT, BUT NOT SURE ITS NEEDED
        this.authKey = UUID.randomUUID().toString();
    }
    
    public String getAuthKey() {
        return this.authKey;
    }
    
    /**
     * When a Person logs out, it gets relayed onto here which will also
     * log out all the Clients.  This is a KISS security solution, but may
     * ultimately need to be more complex.  This will log you out from everywhere.
     */
    public void logOut() {
        for(Client device : clients) {
            device.lockout();
        }
        //Persistor.get().persist(this);
        // new:
        //this.save();
    }
    
    /**
     * If a Person logs in, that gets relayed to the user and then logIn logs in the Client
     * @param crtf 
     */
    public void logIn(Client device) {
        isLoggedIn = true;
        device.login();
        //Persistor.get().persist(this);
        // new:
        //this.save();
    }
    
    private String passwordHash;
    private Set<Client> clients = new HashSet<Client>();
    
    //This is where the login state of the user is stored globally (so you can log yourself out everywhere by setting to false anywhere)
    //HMMM THIS NEEDS MORE THOUGHT, I DON'T KNOW IF THIS IS A PREFERENCE SETTING OR THE ACTUAL STATE
    //THE STATE HAS TO BE HELD FOR THE CLIENT+USER, NOT JUST THE USER AS THEY COULD BE LOGGED IN ON MULTIPLE MACHINES
    //I THINK THIS IS JUST THE GLOBAL SETTING.  WILL NEED TO THINK ABOUT THIS.  MAYBE YOU ARE DOING TWO LOGINS:  ONE FOR A Client
    // and one for the user, and only if both are true are you really logged in.  You need to be able to shut things off locally
    //in case a computer gets taken, so ok, yeah, it's right...this is the state of User login
    private boolean isLoggedIn = false;
    
    //This is where the linkage between the user and its in-clotho Person representations are stored
    private Set<String> personUUIDs = new HashSet<String>();

}
