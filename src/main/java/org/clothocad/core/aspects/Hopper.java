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

import java.util.HashMap;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.ObjectId;

/**
 * The hopper is where Doo's are held while a process is going on elsewhere.
 * 
 * As it's really going to be annoying if Doos are lost due to errors and timeouts,
 * it's important that these things not get dropped.
 * 
 * So, this Hopper thing is for managing this aspect of the Doo's lifecycle.  It needs
 * to be persisted regularly, unlike the Collector where it doesn't matter whether it drops
 * its objects because upon request it just goes and gets them.  Doos are different than
 * Sharables -- they lose their persistance unless you explicitly persist them.  It's in this
 * Hopper state where such losses are likely to occur breaking the chains, so this is the state
 * that needs to be saved.  It should also manage timeouts -- like keep track of how long a Doo
 * has been sitting around waiting, and do something appropriate with it (like logging it, possibly
 * indicating in a log of bad-other-clothos that they don't respond well, stuff like that can be
 * done with a stalled Doo).  Also, you want to time out the stay here or they will accumulate, and
 * since the state of the Hopper needs to be persisted across run-times, it really would fill up with
 * timed-out Doos.
 * 
 * @author John Christopher Anderson
 */

@Slf4j
public class Hopper implements Aspect {    /**
     * Pull the Doo out of the Hopper and return it
     * @param dooId
     * @return 
     */
    public Doo extract(ObjectId dooId) {
        Doo extracted = dooList.get(dooId);
        dooList.remove(dooId);
        return extracted;
    }
    
    public Doo extract(String dooId) {
        return extract(new ObjectId(dooId));
    }
    
    /**
     * For an operation that involves putting a Doo in cold-storage, when later
     * you revive that Doo, you have two options -- you can terminate its
     * lifecycle (with this method), or you can just pop it out (extract).
     * 
     * @param dooId
     * @return 
     */
    public boolean terminate(ObjectId dooId) {
        Doo doo = dooList.get(dooId);

        if(doo==null) {
            //THIS SHOULD BE AN ERROR-MESSAGE SCENARIO FOR SURE AS IT DROPPED A DOO, WHICH IS A NO-NO FOR SURE.
            log.error("terminate: Cannot find doo with id {}", dooId);
            return false;
        }
        //IF THE DOO TELLS CLOTHO TO DO ANYTHING NOW, DOO IT HERE
        
        //Terminate this process
        doo.setMessage("Process successfully completed.");
        doo.terminate();
        dooList.remove(dooId);
        return true;
    }
    
    /**
     * Add a Doo to the Hopper
     * @param doo 
     */
    public void add(Doo doo) {
        dooList.put(doo.getId(), doo);
    }
    
    //Singleton stuff
    private Hopper() { }
    private static final Hopper singleton = new Hopper();
    public static Hopper get() {
        return singleton;
    }
    
    
    private HashMap<ObjectId, Doo> dooList = new HashMap<>();
    
    
}
