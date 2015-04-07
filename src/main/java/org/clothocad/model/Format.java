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

package org.clothocad.model;

import com.fasterxml.jackson.annotation.JsonTypeInfo;
import java.util.List;

/**
 * Interface for a ClothoFormat
 * @author J.Christopher Anderson
 * @author Douglas Densmore
 * @author Nicholas Roehner
 */ 
@JsonTypeInfo(use = JsonTypeInfo.Id.CLASS, property = "schema", include = JsonTypeInfo.As.PROPERTY)
public interface Format {

    /**
     * Test whether a part is a valid member of the ClothoFormat specs
     * @param p Part to test
     * @return true if adheres to the format; false otherwise
     */
    public boolean checkPart(Part p);

    /**
     * Check to see if an arraylist of parts permits making a composite part
     * @param composition
     * @return
     */
    public boolean checkComposite(List<Part> subParts);

     /**
     * Generates composite part in accordance with this format
     * @param name
     * @param composition
     * @param author
     * @return
     */
    public Part generateCompositePart(String name, List<Part> subParts, Person author);
    
    /**
     * Generates composite part in accordance with this format
     * @param name
     * @param description
     * @param composition
     * @param author
     * @return
     */
    public Part generateCompositePart(String name, String description, List<Part> subParts, Person author);
  
}
