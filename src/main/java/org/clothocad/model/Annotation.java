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

import org.clothocad.core.datums.ObjBase;

import lombok.Getter;
import lombok.NoArgsConstructor;
import org.clothocad.core.persistence.annotations.Reference;

/**
 * An Annotation is a single line of genbank essentially.  It maps a Feature
 * or just something a user has labeled as signficant to a NucSeq.
 */

@NoArgsConstructor
public class Annotation extends ObjBase {

    /**
     * Constructor for an Annotation that is not a Feature, just a region of colored sequence
     *
     * @param name
     * @param nucseqid  the NucSeq that you're annotating
     * @param forColor
     * @param revColor
     * @param Start
     * @param End
     * @param user
     * @param plusstrandtrue
     * @param symbol  (can be null)
     */
    public Annotation( String name, NucSeq nucseq, Color forColor, Color revColor, int start, int end, Person user, boolean plusstrandtrue, String symbol ) {
        super(name);

        sequence = nucseq;
        forwardColor = forColor;
        reverseColor = revColor;
        this.start = start;
        this.end = end;
        author = user;
        isForwardStrand = plusstrandtrue;
        this.symbol = symbol;
        nucseq.addAnnotation(this);
    }

    /**
     * Constructor for an Annotation that corresponds to a Feature object
     * @param afeature
     * @param nucseqid
     * @param forColor
     * @param revColor
     * @param Start
     * @param End
     * @param user
     * @param plusstrandtrue
     * @param symbol
     */
    public Annotation( Feature afeature, NucSeq nucseq, Color forColor, Color revColor, int start, int end, Person user, boolean plusstrandtrue, String symbol ) {
        this(afeature.getName(), nucseq, forColor, revColor, start, end, user, plusstrandtrue, symbol);
        feature = afeature;
        if ( forColor == null ) {
            forwardColor = afeature.getForwardColor();
        }
        else {
            forwardColor = forColor;
        }
        if ( revColor == null ) {
            reverseColor = afeature.getReverseColor();
        }
        else{
            reverseColor = revColor;
        }
    }



    /**
     * Reverse the orientation of the annotation (reverse complement
     * it and flip flop the start and end sites).  Called from NucSeq
     * when it's reverse complemented
     * @param nucseqLength
     */
    void invert(int nucseqLength){
        isForwardStrand = !isForwardStrand;
        int s = start;
        start = nucseqLength - end;
        end = nucseqLength - s;
    }

    /**
     * Get the approriate color for the annoation
     * @return either the forward or reverse color depending
     * on the orientation of the annotation
     */
    public Color getColor() {
        if ( isForwardStrand ) {
            return forwardColor;
        }
        return reverseColor;
    }

    /**
     * Get the forward color as an integer code
     * @return an integer of the Color
     */
    public int getForwardColorAsInt() {
        return forwardColor.getRGB();
    }
    /**
     * Get the reverse color as an integer code
     * @return an integer of the Color
     */
    public int getReverseColorAsInt() {
        return reverseColor.getRGB();
    }

    /*-----------------
    variables
    -----------------*/
    @Getter
    private String symbol;

    @Getter
    private boolean isForwardStrand;
    @Getter
    @Reference
    private Person author;
    @Getter
    @Reference
    private Feature feature;
    @Getter
    @Reference
    private NucSeq sequence;
    @Getter
    private int start, end;
    @Getter
    private Color forwardColor, reverseColor;
    

}
