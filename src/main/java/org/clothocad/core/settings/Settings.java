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

package org.clothocad.core.settings;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.Sharable;

/**
 * The settings wraps a Datum to store the settings.
 * 
 * The philosophy here is to fill this up globally with whatever is
 * convenient for coding things that need to access this info.  That
 * will make for highest retrieval rate during gets which will be far
 * more common than sets. The setters can be written later as the
 * pattern of things that would make sense to a user such that it need
 * not be so verbose.  So, intentions will be translated into these
 * more specific things.
 * 
 * Doing it as a persisted class will guarantee default settings are
 * entirely conservative on all fronts.
 * @author John Christopher Anderson
 */
public class Settings {
    public static boolean isRecordAllDoos() {
        return datum.recordAllDoos;
    }
    
    public static String getRootURL() {
        //THIS SHOULD BE REPLACED BY A DYNAMIC DETERMINATION OF THE CURRENT IP
        return "http://" + getHost() +
               ((getPort()==80) ? "" : ":" + getPort().toString());
    }

    public static String getHost() {
        return "localhost";
    }

    public static Integer getPort() {
        return 8080;
        //return 61623;
    }

    private class SettingsDatum 
    		extends ObjBase {
        private boolean recordAllDoos = false;
        private String id = "settings-datum-is-uuid";
    }
    
    private static SettingsDatum datum;
}
