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

package org.clothocad.core.communication;

import java.lang.ref.WeakReference;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.communication.mind.Widget;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.util.JSON;

/**
 * @author John Christopher Anderson
 */


public class UpdateIndex {
    Persistor persistor;

    public void register(Sharable sharable, Widget widget, Mind mind) {
        //Register the sharable --> widget hash
        ObjectId sharableId = sharable.getId();
        Set<Widget> widgets = null;
        if(sharableToWidget.containsKey(sharableId)) {
            widgets = sharableToWidget.get(sharableId);
        } else {
            widgets = new HashSet<>();
        }
        widgets.add(widget);
        sharableToWidget.put(sharableId, widgets);
        
        //Register the widget --> mind hash
        widgetToMind.put(widget.getId(), new WeakReference(mind));
    }
    
    public void update(String sharableId) {
//        Sharable sharable = persistor.get(Sharable.class, new ObjectId(sharableId));
//        update(sharable);
    }
    
    public void update(Sharable sharable) {
        Set<Widget> widgets = sharableToWidget.get(sharable.getId());
        if(widgets==null) {
            System.out.println("UpdateIndex.update(...) has null for widgets, so aborting.");
            return;
        }
        
        Map<String, Object> data = JSON.mappify(sharable);
        
        for(Widget widg : widgets) {
            try {
                String socket_id = widg.getPage().getSocketId();
                Mind mind = widgetToMind.get(widg.getId()).get();
            } catch(Exception err) {
                err.printStackTrace();
            }
        }
    }

    private HashMap<ObjectId, Set<Widget>> sharableToWidget = new HashMap<>();
    private HashMap<String, WeakReference<Mind>> widgetToMind = new HashMap<String, WeakReference<Mind>>();
}
