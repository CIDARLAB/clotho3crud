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
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.aspects;

import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.clothocad.core.communication.Callback;

/**
 * @author John Christopher Anderson
 */

public class Process {
    
    public static void nextTick(Callback callback, Exception err) {
        singleton.executeOnError(callback, err);
    }
    
    public static void nextTick(Callback callback, Map<String, Object> data) {
        singleton.executeOnSuccess(callback, data);
    }
    
    private void executeOnSuccess(final Callback callback, final Map<String, Object> data) {
        executorService.submit(new Runnable(){
            @Override
            public void run() {
                callback.onSuccess(data);
            }
        });
    }
    
    private void executeOnError(final Callback callback, final Exception err) {
        executorService.submit(new Runnable(){
            @Override
            public void run() {
                callback.onFailure(err);
            }
        });
    }

    public static Process get() {
        return singleton;
    }
    
    public final ExecutorService executorService = Executors.newFixedThreadPool(1);
    private static final Process singleton = new Process();
}
