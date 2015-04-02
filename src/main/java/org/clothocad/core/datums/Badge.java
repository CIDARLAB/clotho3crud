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

import java.util.Map;
import javax.validation.constraints.AssertTrue;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.model.Person;


/**
 * Badges are visually apparent as animated gifs supplied by Trails to indicate
 * training or as prizes for the sharing of large amounts of Data.
 * 
 * So, they are used for motivation, biosafety, and so people can better understand
 * the role of individual Persons.
 * 
 * Badges live as instances in the container.  They aren't defined by a Schema because
 * of security reasons they need to stay out of the container.
 * 
 * Only the Core can assign a Badge to a Person.
 * 
 * They definitely aren't Sharables but are intrinsically always publically viewable
 * and creatable from within the container (for which they'll need their own editor)
 * 
 * The fields here should be DataFields so that this can be indexed by normal semantic types
 * (This isn't what's implemented at present)
 * 
 * @author John Christopher Anderson
 */


public class Badge 
		extends SharableObjBase{

	public Badge(Person author) {
		super("", author);
	}

	//
//    @Override
//    public boolean canView(String userId) {
//        //Badges are always viewable by anybody
//        return true;
//    }
//
//    @Override
//    public boolean canEdit(String userId) {
//        //The creator of the Trail can edit it
//        if(userId.equals(creator.get().getUUID())) {
//            return true;
//        }
//        
//        //If they are in my list of editors they can edit it
//        if(this.authorizedEditors.contains(userId)) {
//            return true;
//        }
//        
//        //Otherwise they can't
//        return false;
//    }
//
//    @Override
//    public boolean canShare(String userId) {
//        //State maintenance of a Badge is never sent elsewhere
//        return false;
//    }
//
//    @Override
//    public boolean canSell(String userId) {
//        return false;
//    }
//
//    @Override
//    public Instance getOwner() {
//        return (Instance) this.creator.get();
//    }
//
//    @Override
//    public final Map<String, Object> getJSON() {
//        try {
//            Map<String, Object> out = new Map<String, Object>();
//            out.put("type", "BADGE");
//            out.put("name", this.name);
//            out.put("description", description);
//            out.put("smallIcon", (String) this.smallIcon.getValue());
//            out.put("largeIcon", (String) this.largeIcon.getValue());
//            out.put("authURL", Settings.getRootURL());
//            return out;
//        } catch(Exception err) {
//            return null;
//        }
//    }
//
//    @Override
//    public void set(Map<String, Object> permissions, String userId) {
//        
//    }
//
//    @Override
//    public SharableType getSharableType() {
//        return SharableType.BADGE;
//    }
//    
    public final boolean hasBadge(String personId) {
return true;
        //        return haveBadge.contains(personId);
    }

	//XXX: move to ObjBase or Sharable
    @AssertTrue
    public boolean validate(Map<String, Object> obj) {
        return true;
    }
//    
    /**
     * Will need some stuff in here to denote popularity and uniqueness...it needs to keep track somehow
     * of how meaningful these badges are so it can sort them appropriately and automatically highlight
     * individuals one day
     */
    
    
    /**
     * Hmmm, in terms of functionality, somewhere you need to define the badge's class
     * as being the 'golden pipette badge' meaning you've shared over 1000 parts or
     * or the equivalent.
     * 
     * Errr, no.  It's there are these Badge objects that have THIS form, and people get on their list...yeah.
     * 
     * These need to be editable from within the container but aren't sharable and such (they are intrinsically viewable but never sharable, and editing is restricted to a list of users)
     */
    
//    private IconField smallIcon; //The animated gif in a format convenient for web
//    private IconField largeIcon; //A larger (or more detailed) version of the icon
//    private NameField name;
//    private StringField description;
//    
//    private XRefBag haveBadge; //The list of uuids of Persons that have earned this badge
//    private AuthorField creator;
//    private XRefBag authorizedEditors;
}
