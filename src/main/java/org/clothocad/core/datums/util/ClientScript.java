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

package org.clothocad.core.datums.util;

import org.clothocad.core.datums.ObjBase;

/* This one deals with passing json:
 * http://stackoverflow.com/questions/1078419/java-scriptengine-using-value-on-java-side
 * http://docs.oracle.com/javase/6/docs/technotes/guides/scripting/programmer_guide/index.html
 */

/**
 * @author John Christopher Anderson
 */
public class ClientScript extends ObjBase {
    private ClientScript() { }
    
    public ClientScript(String js, String html, Language language) {
        this.onShowScript = js;
        this.graphicsCode = html;
        this.language = language;
    }

    public Language getLanguage() {
        return language;
    }

    public String getGraphicsCode() {
        return graphicsCode;
    }

    public String getOnShowScript() {
        return onShowScript;
    }
    
    private String onShowScript;
    private String graphicsCode;
    private Language language;
}
