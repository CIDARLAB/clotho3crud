/*
Copyright (c) 2009 The Regents of the University of California.
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

import java.awt.Color;

import javax.validation.constraints.NotNull;

import org.clothocad.core.datums.SharableObjBase;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.persistence.annotations.Reference;

/**
 * An Annotation is a single line of Genbank essentially.  It maps a Feature
 * or just something a user has labeled as signficant on a Sequence.
 */

/**
 * 
 * @author J. Christopher Anderson
 * @author Nicholas Roehner
 */

public class Annotation extends SharableObjBase {

	@Getter
	@Setter
	protected String symbol;
    
	@NotNull
	@Getter
    @Setter
    protected boolean isForwardStrand;
    
    @Getter
    @Setter
    @Reference
    protected Feature feature;
    
    @NotNull
    @Getter
    @Setter
    protected int start, end;
    
    @Setter
    protected Color forwardColor, reverseColor;
	
    /**
     * @param name
     * @param seq
     * @param end
     * @param start
     * @param author
     * @param isForwardStrand
     */
   protected Annotation(String name, int start, int end, boolean isForwardStrand, Person author) {
        super(name, author);
        this.start = start;
        this.end = end;
        this.isForwardStrand = isForwardStrand;
    }
    
    /**
     * @param name
     * @param description
     * @param seq
     * @param end
     * @param start
     * @param author
     * @param isForwardStrand
     */
   protected Annotation(String name, String description, int start, int end, 
    		boolean isForwardStrand, Person author) {
        super(name, author, description);
        this.start = start;
        this.end = end;
        this.isForwardStrand = isForwardStrand;
    }

    /**
     * Reverse the orientation of the annotation (reverse complement
     * it and flip flop the start and end sites).  Called from NucSeq
     * when it's reverse complemented
     * @param seqLength
     */
    public void invert(int seqLength){
        isForwardStrand = !isForwardStrand;
        int s = start;
        start = seqLength - end;
        end = seqLength - s;
    }

    /**
     * Get the approriate color for the annoation
     * @return either the forward or reverse color depending
     * on the orientation of the annotation
     */
    public Color getColor() {
        if (isForwardStrand) {
            return getForwardColor();
        } else {
        	return getReverseColor();
        }
    }

    /**
     * Get the forward color as an integer code
     * @return an integer of the Color
     */
    public int getForwardColorAsInt() {
        return getForwardColor().getRGB();
    }
    /**
     * Get the reverse color as an integer code
     * @return an integer of the Color
     */
    public int getReverseColorAsInt() {
        return getReverseColor().getRGB();
    }
    
    /**
     * Get the preferred forward color for this Annotation.  If no forward color
     * was set, a default color will be returned.
     * @return an AWT Color object.  It won't be null;
     */
    public Color getForwardColor() {
        if (forwardColor == null) {
        	forwardColor = new Color(125, 225, 235);
        }
        return forwardColor;
    }

    /**
     * Get the preferred reverse color for this Annotation.  If no reverse color
     * was set, a default color will be returned.
     * @return an AWT Color object.  It won't be null;
     */
    public Color getReverseColor() {
        if (reverseColor == null) {
        	reverseColor = new Color(125, 225, 235);
        }
        return reverseColor;
    }
    
    /**
     * Set the forward and reverse preferred colors for this feature to some
     * random medium-bright color.
     */
    public void setRandomColors() {
    	int[][] intVal = new int[2][3];
    	for (int j = 0; j < 3; j++) {
    		double doubVal = Math.floor(Math.random() * 155 + 100);
    		intVal[0][j] = (int) doubVal;
    		intVal[1][j] = 255 - intVal[0][j];
    	}
    	forwardColor = new Color(intVal[0][0], intVal[0][1], intVal[0][2]);
    	reverseColor = new Color(intVal[1][0], intVal[1][1], intVal[1][2]);
    }

}
