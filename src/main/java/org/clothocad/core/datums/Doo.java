/*
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
ENHANCEMENTS, OR MODIFICATIONS.
 */

package org.clothocad.core.datums;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import org.clothocad.core.settings.Settings;

/**
 * @author John Christopher Anderson
 */
public class Doo 
		extends ObjBase {
    public Doo(Doo parent, boolean saveit) {
        if(parent == null) {
            parentDooId = null;
        } else {
            parentDooId = parent.getId();
        }
        savePolicy = saveit;
        save();
    }
    
    public void setMessage(String str) {
        MsgLine msg = new MsgLine();
        msg.msg = str;
        messages.add(msg);
        save();
    }
    
    public String getMessage() {
        if(messages.isEmpty()) {
            return "empty";
        }
        return messages.get(messages.size()-1).msg;
    }

    public void abort(Exception err) {
        setMessage("Doo aborted due to Exception");
        abortErr = err;
        terminate();
    }
    
    public void terminate() {
        setMessage("Doo terminated by Doo.terminate()");
        dateEnded = new Date();
      //if global settings say to....
        //  Persistor.get().save(this);
    }
    
    public boolean save() {
        if(!savePolicy) {
            return false;
        }
        
        if(Settings.isRecordAllDoos()) {
            //Persistor.get().persistDatum(this);
        	this.save(); //why infinite loop?
        }
        return true;
    }
    
    private static class MsgLine {
        String msg;
        Date date;
    }
    
    private final ObjectId parentDooId;    
    private List<MsgLine> messages = new ArrayList<MsgLine>();
    private Date dateEnded = null;
    private boolean savePolicy = false;
    private Exception abortErr;

}
