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

package org.clothocad.core.aspects.Interpreter;

import java.util.Map;
import java.util.TreeMap;
import java.util.Set;

/**
 * @author John Christopher Anderson
 */
public class Interpreter {
    /**
     * Entry point B.
     * Receives nlp commands such as 'run pBca1256 on sequenceview'
     * Ultimately this either needs a generic callback or a session reference
     * or something so this can all respond.
     * 
     * @param cmd 
     */
    public Set<String> receiveNative(String cmd) {
        String[] features = Handler.convertToFeatures(Utilities.tokenize(cmd));
        TreeMap<String, Double> SortedActionList = Handler.query(features);
        return SortedActionList.keySet();
    }

    /**
     * Entry point that communicates with learner for learning 
     * the association of a command with its action statement.
     * @param cmd
     * @param action
     */
    public void learnNative(String cmd, Map<String, Object> action) {
        Learner.learn(cmd, action);
    }
    
    public void forgetNative(String cmd, String action) {
        Learner.forget(cmd, action);
    }

    public static Interpreter get() {
        return singleton;
    }

    private static final Interpreter singleton = new Interpreter();
}
