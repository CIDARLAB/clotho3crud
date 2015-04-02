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

import org.clothocad.core.persistence.annotations.Reference;
import java.awt.Color;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.SharableObjBase;

/**
 *
 * @author J. Christopher Anderson
 */

//unique name requirement
//sequence should be unique
//start/stop codons get clipped
//sequence cannot be degenerate
//cds-like sequences might be cds
//cds features must have valid cds sequence

//took out notes for demo

@NoArgsConstructor
public class Feature extends SharableObjBase {
    /**
     * Relayed constructor of a new Feature
     * @param name
     * @param seq
     * @param author
     */
    private Feature(String name, NucSeq seq, Person author, boolean iscds) {
        super(name, author);
        sequence = seq;
        isCDS = iscds;
    }

    /**
     * Method for generating a Feature.
     * Features are either CDS's or not.  For non-CDS features, the Feature will have the sequence as
     * explicitly submitted to the method.  For CDS features, you can provide a sequence with start and
     * stop codons on the ends, but they'll get chopped off.  You can also provide them without starts or stops
     * as long as they are in frame from start to end.  Starts and stops are inferred later during
     * autoannotate, so there is no real information loss.  The reason things are this way is to avoid calling
     * multiple subtly-different CDS's different names when they are simply regulatory variants of the same
     * polypeptide-coding sequence.  Different codon usage, however, is regarded as being a separate feature
     * as it may behave differently during translation (folding rates, translation rates, etc.)
     *
     * Regardless of strand, a Feature is considered to be equivalent.  The CDS Features must be passed as if
     * they are encoded on the sense strand (ATG...TAA).  All Features will read 5'
     * to 3', though they can be annotated onto NucSeqs in either orientation.
     *
     * For CDS parts with intentional internal stops (often used during amber suppression, but also might
     * be natural for some organisms), this method will warn the user that they are present, but it is acceptable
     * to have them.
     *
     * If a sequence passed to this method has degeneracy codes (N's R's etc.), it currently will be rejected.
     * Future version of Clotho will support libraries, but for now, Features are exact sequences.
     *
     * After creating a Feature, you can set the Forward and Reverse colors programmatically.  If none are
     * supplied, random colors will be generated when those colors are later requested.
     *
     * You can also add Notes to the Feature after creating the Feature.
     *
     * Setting the Family(s) for a Feature is currently not implemented, but will be implemented in future versions
     * of Clotho.
     *
     * @param name  the Feature's name (should be unique in the database)
     * @param seq the sequence of the Feature (case is ignored)
     * @param author  the Person object to be author of the new Feature
     * @param type
     * @return the new Feature object, a preexisting Feature of the same sequence, or null if the new Feature
     * was rejected for some reason.
     */
    public static Feature generateFeature(String name, String seq, Person author, boolean iscds) {
        // To find a feature who's sequence matches the above (chop off start and stop for CDS)
        String uppSeq = seq.toUpperCase();
        String testseq = uppSeq;
        if (uppSeq.startsWith("ATG") || uppSeq.startsWith("GTG")) {
            testseq = uppSeq.substring(3);
        }
        if (uppSeq.endsWith("TAA") || uppSeq.endsWith("TAG") || uppSeq.endsWith("TGA")) {
            String dudly = testseq.substring(0, testseq.length() - 3);
            testseq = dudly;
        }

        //If it's cleared checks for preexisting features, start creating a new one
        NucSeq nseq = new NucSeq(uppSeq);

       /* //If it starts with a start codon, but isn't stated as a CDS, maybe throw a warning
        if (!iscds) {
            //Check if it's a CDS in forward orientation
            if (uppSeq.startsWith("ATG")) {
                int modulus = nseq.getSeq().length() % 3;
                int extras = nseq.getSeq().length() - modulus;
                String translation = nseq.translate(0, extras);
                if (translation.equals("*")) {
                    System.out.println("Failed translation of " + nseq.getSeq().substring(0, extras));
                }
                //If the first stop codon is near the end, it might be a CDS
                if (translation.indexOf("*") > translation.length() - 3 || translation.indexOf("*") == -1) {
                    int n = JOptionPane.showConfirmDialog(
                            null,
                            "Feature " + name + " might be a CDS Feature, is it?",
                            "CDS detected",
                            JOptionPane.YES_NO_OPTION);
                    if (n == 0) {
                        iscds = true;
                    }
                }
            }
            //Check if it's a CDS in reverse orientation
            String srevcomp = nseq.revComp();
            NucSeq revcomp = new NucSeq(srevcomp);
            revcomp.setTransient();
            if (srevcomp.startsWith("ATG")) {
                int modulus = revcomp.getSeq().length() % 3;
                int extras = revcomp.getSeq().length() - modulus;
                String translation = revcomp.translate(0, extras);
                if (translation.equals("*")) {
                    System.out.println("Failed translation of " + nseq.getSeq().substring(0, extras));
                }
                //If the first stop codon is near the end, it might be a CDS
                if (translation.indexOf("*") > translation.length() - 3 || translation.indexOf("*") == -1) {
                    int n = JOptionPane.showConfirmDialog(
                            null,
                            "Feature " + name + " might be a CDS Feature reverse complemented, is it?",
                            "CDS detected",
                            JOptionPane.YES_NO_OPTION);
                    if (n == 0) {
                        iscds = true;
                        nseq.changeSeq(revcomp.getSeq());
                    }
                }
            }
        }*/

        /*
        //If it was stated to be a CDS, make sure it's really a CDS
        if (iscds) {
            System.out.println("doing iscds logic checks");
            String translation = nseq.translate(0);
            //If there's a start codon on it, chop it off
            if (uppSeq.startsWith("ATG") || uppSeq.startsWith("GTG") || uppSeq.startsWith("TTG")) {
                System.out.println("I'm chopping the start codon");
                String tempseq = uppSeq.substring(3);
                uppSeq = tempseq;
                nseq.changeSeq(uppSeq);
                System.out.println("done chopping the start codon");
            }

            //Trim the end if it's out of frame
            int extras = uppSeq.length() % 3;
            if (extras > 0) {
                nseq.changeSeq(uppSeq.substring(0, uppSeq.length() - extras));
            }

            //If there's a stop codon at the end, chop it off
            translation = nseq.translate(0);
            int nearend = translation.length() - 2;
            int end = translation.indexOf("*");
            System.out.println(translation + " near end " + nearend + " end " + end);
            if (end > nearend) {
                System.out.println("chopping off the stop codon");
                nseq.changeSeq(uppSeq.substring(0, 3 * (end)));

                //Otherwise it doesn't translate through, so maybe dump it
            } else if (end == -1) {
                System.out.println("there was no stop codon");
            } else {
                int n = JOptionPane.showConfirmDialog(null, "Feature " + name + " appears to contain an internal start codon.  Should I cancel?", "Internal stop "
                        + "detected", JOptionPane.YES_NO_OPTION);
                if (n == 0) {
                    nseq.setTransient();
                    return null;
                }
            }

            //If the thing translates wrong, return null
            if (nseq.translate(0).equals("*")) {
                nseq.setTransient();
                return null;
            }
        }
        */

        final Feature afeature = new Feature(name, nseq, author, iscds);

       /* //Set the biosafety level of the new Feature
        Thread bslThread = new Thread() {

            @Override
            public void run() {
                afeature.setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.NAME_CHANGED);
                Short rg = afeature.getSeq().performBiosafetyCheck();
                afeature._featDatum._riskGroup = rg;
            }
        };
        bslThread.start();
        afeature.addSaveHold(bslThread);*/

        //TODO: validate
        
        nseq.setLocked(true);
        //   NucSeq.addFeatureToTable(afeature);
        return afeature;
    }

    /**
     * This method may exist in future versions of Clotho but is not currently implemented.
     *
     * Method for altering a Feature that is a mutant of the first.  The most likely scenario for this is
     * when a Feature already exists and is linked to a part, but a mutant plasmid is generated, and the author
     * needs to indicate the existence of the mutant and instead of calling it a new Feature, just allow for it
     * to be a little degenerate.
     * 
     * @param name
     * @param seq
     * @param author
     * @return
     *
    public boolean mutateFeature( String name, String seq, Person author ) {
    String old = this.getSeq().toString();
    String nu = "";
    int Ncount = 0;
    
    for ( int i = 0; i < old.length(); i++ ) {
    if ( old.charAt( i ) == seq.charAt( i ) ) {
    nu += old.charAt( i );
    } else {
    nu += 'N';
    Ncount++;
    }
    }
    if ( Ncount > 25 ) {
    JOptionPane.showMessageDialog( null, "That sequence is pretty far from the original, I count " + Ncount + " differences.  You should make a new feature instead, or check the alignment of the new sequence.", "Error", JOptionPane.ERROR_MESSAGE );
    return false;
    }
    _featDatum._seqID = new NucSeq( nu ).getId();
    return true;
    }
     */

    /**
     * Get the preferred forward color for this Feature.  If no forward color
     * was set, a default color will be returned.
     * @return an AWT Color object.  It won't be null;
     */
    public Color getForwardColor() {
        if (forwardColor == null) {
            return new Color(125, 225, 235);
        }
        return forwardColor;
    }

    /**
     * Get the preferred reverse color for this Feature.  If no reverse color
     * was set, a default color will be returned.
     * @return an AWT Color object.  It won't be null;
     */
    public Color getReverseColor() {
        if (reverseColor == null) {
            return new Color(125, 225, 235);
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

    /**
     * Retrieve a Feature from the database by its name
     * @param name the name of the desired Feature
     * @return the Feature or null if none was found
     */
    public static Feature retrieveByName(String name) {
        throw new UnsupportedOperationException();
    }


    /**
     * General method for adding things to the Feature.  Currently it only
     * accepts additions of Notes.
     * @param dropObject the ObjBase being added (linked) to this Feature
     * @return true if the drop was accepted
     */
    /*@Override
    public boolean addObject(ObjBase dropObject) {
        switch (dropObject.getType()) {
            case NOTE:
                final Note item = (Note) dropObject;

                ActionListener undo = new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        _featDatum._noteLinks.remove(item.getId());
                        item.removeFeatureRelay(Feature.this.getId());
                        setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.NOTE_LINKED);
                    }
                };
                ActionListener redo = new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        _featDatum._noteLinks.add(item.getId());
                        item.addFeatureRelay(Feature.this);
                        setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.NOTE_LINKED);
                    }
                };
                addUndo(undo, redo);

                _featDatum._noteLinks.add(item.getId());
                item.addFeatureRelay(this);
                setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.NOTE_LINKED);
                return true;
            case FAMILY:
                final Family fitem = (Family) dropObject;

                System.out.println("Addding Family " + fitem.getName() + " to " + Feature.this.getName());
                ActionListener fundo = new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        _featDatum._familyLinks.remove(fitem.getId());
                        fitem.removeFeatureRelay(Feature.this.getId());
                        setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.FAMILY_TO_FEATURE);
                    }
                };
                ActionListener fredo = new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent e) {
                        _featDatum._familyLinks.add(fitem.getId());
                        fitem.addFeatureRelay(Feature.this);
                        setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.FAMILY_TO_FEATURE);
                    }
                };
                addUndo(fundo, fredo);

                _featDatum._familyLinks.add(fitem.getId());
                fitem.addFeatureRelay(this);
                setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.FAMILY_TO_FEATURE);
                return true;
            default:
                return false;
        }
    }*/


    /* SETTERS
     * */
    /**
     * Method for changing the sequence of this Feature.  The sequence of the NucSeq is
     * locked in NucSeq.  You must use this method to alter the sequence instead.
     * @param newseq the new sequence of the Feature
     */
    public void setSequenceSequence(final String newseq) {
        if (newseq == null || newseq.equals("")) {
            //fireData(new RefreshEvent(this, RefreshEvent.Condition.SEQUENCE_CHANGED));
            return;
        }
        sequence.APIchangeSeq(newseq);

        //todo: Change the risk group
    }

    /**
     * Change the risk group of the Feature.  You can only raise the risk group.
     * @param newrg the new risk group (1 through 5)
     */
    public final void setRiskGroup(short newrg) {
        if (newrg > riskGroup && newrg <= 5) {
            //addUndo("_riskGroup", _featDatum._riskGroup, newrg);
            riskGroup = newrg;
           // setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.RISK_GROUP_CHANGED);
        }
        //todo: throw appropriate invalid operation exception
    }

    /* GETTERS
     * */
    /**
     * Get the Author of this Feature as a UUID link.
     * @return a UUID String
     */
    public ObjectId getAuthorUUID() {
        return author.getId();
    }



    /*-----------------
    variables
    -----------------*/
    @Setter
    @Getter
    @Reference
    private NucSeq sequence;
    
    @Setter
    private Color forwardColor, reverseColor;
    @Setter
    @Getter
    private String genbankId, swissProtId, PDBId;
    //private String sourceOrganism; was in datum
    
    @Setter
    @Getter
    @Reference
    private Person author;
    @Getter
    private short riskGroup;
    //private String featureData; was in datum

       
    @Getter
    private boolean isCDS;

}
