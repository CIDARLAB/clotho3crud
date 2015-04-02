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
createQuery
THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */
package org.clothocad.model;

import org.clothocad.core.persistence.annotations.Reference;
import java.util.List;
import java.util.Map;
import javax.validation.constraints.AssertTrue;
import javax.validation.constraints.NotNull;

import org.clothocad.core.datums.ObjBase;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import lombok.ToString;
import lombok.extern.slf4j.Slf4j;

/**
 *
 * @author jcanderson
 */
@ToString(callSuper=true, includeFieldNames=true)
@NoArgsConstructor
@Slf4j
public abstract class Part extends ObjBase {

    public static Part generateBasic(String name, String description, String seq, Format format, Person author) {
        return new BasicPart(name, description, seq, format, author);
    }
    
    //inject connection

    public static Part generateComposite(List<Part> composition, Object additionalRequirements, Format format, Person author, String name, String description) {
        return new CompositePart(composition, additionalRequirements, format, author, name, description);
    }

    
    @Setter
    @Getter
    @Reference
    private Person author;
    
    @Getter
    @NotNull
//XXX:    @Replace(encoder = "getFormatName", decoder="setFormatFromName")
    private Format format;

    @Getter
    @Setter
    private String shortDescription;
    
    @Getter
    @Setter
    private PartFunction type;
    
    @Getter
    private short riskGroup;
    
    protected Part(String name, String desc, Format form, Person author){
        super(name); 
        this.shortDescription = desc;
        this.author = author;
        this.format = form;

    }

    /**
     * Create composite from scratch
     *
     * @param name  nickname of Part, like "roo40"
     * @param shortdescription short description, such as "[TetR]"
     * @param seq sequence of the Part like "cgaaggcaggacacacatg"
     * @param form Part Format
     * @param author author of the Part
     * @param partType Basic or Composite
     */
    /**
     * Call this method to construct a new basic Part.  It will check that:
     *  1) A sequence in this Format isn't already in the database
     *      otherwise it returns the Part that already exists
     *  2) The Part obeys its Format standard
     *      otherwise it returns null
     *
     * Only if those are satisfied is it added to the Collector and returned to
     * the calling code.
     *
     * @param name  nickname of Part, like "roo40"
     * @param shortdescription short description, such as "[TetR]"
     * @param seq sequence of the Part like "cgaaggcaggacacacatg"
     * @param form Part Format
     * @param author author of the Part
     * @param partType Basic or Composite
     */
   
    @AssertTrue
    public abstract boolean checkFormat();
    
    public abstract PartType getPartType();
    
    public String getFormatName(){
        return format.getClass().getSimpleName();
    }

    public void setFormatFromName(Map<String, Object> dbObject){
        String name = (String) dbObject.get("format");
        try {
            //XXX: stupid hack
            //TODO: search schemas instead
            format = (Format) Class.forName("org.clothocad.model." + name).newInstance();
            
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException ex) {
            log.error("Couldn't find format {}", "org.clothocad.model." + name);
        } 
    }
    
    /* SETTERS
     * */



    /**
     * Change the Format of the Part
     * @param f
     */
    public void setFormat(Format f) {
        if (f == null) {
            return;
        }

        boolean ok = f.checkPart(this);
        if (!ok) {
            return;
        }

        format = f;
    }


    /* GETTERS
     * */
    /*public List<ObjBase> getPlasmids() {
        ClothoConnection c = Collector.getDefaultConnection();
        ClothoQuery mainQuery = c.createQuery(ObjType.PLASMID);
        ClothoQuery partQuery = mainQuery.createAssociationQuery(Plasmid.Fields.PART);
        partQuery.add(partQuery.getMatchesCrit(Part.Fields.NAME, this.getName()));
        List<ObjBase> results = mainQuery.getResults();
        return results;
    }*/


    public abstract NucSeq getSequence();

    /*public final void changeRiskGroup(Short newrg) {
        if (newrg > _partDatum._riskGroup) {
            addUndo("_riskGroup", _partDatum._riskGroup, newrg);
            _partDatum._riskGroup = newrg;
        }
        setChanged(RefreshEvent.Condition.RISK_GROUP_CHANGED);
    }

    /*
     * Determines the risk group of the Part,
     * relayed from the initial call to NucSeq's
     * call to foreign server
     * @param rg
     *
    private void setRiskGroup(short rg) {
        if (rg == 5) {
            relayRiskGroup((short) 5);
            return;
        }
        if (this._partDatum._partType.equals(partType.Composite)) {
            short currentHighest = rg;
            boolean firsthigher = false;
            for (Part p : this.getCompositeParts()) {
                //If the Part's risk group hasn't been determined, this one isn't either
                if (p.getRiskGroup() == -1) {
                    relayRiskGroup((short) -1);
                    return;
                }

                //If a subpart has a 2+ risk group, increment highest
                if (p.getRiskGroup() > currentHighest) {
                    currentHighest = p.getRiskGroup();
                    relayRiskGroup(currentHighest);
                }

                
                //If a subpart has a 2+ risk group
                if (p.getRiskGroup() > 1) {
                    if (firsthigher) {
                        
                        //Throw a dialog asking for user to put in the new risk group
                        ButtonGroup group = new javax.swing.ButtonGroup();
                        String msgString = "This composite part joins two subparts with risk groups of 2 or higher.  What should the new value be?";
                        int numelements = 5 - currentHighest;
                        Object[] array = new Object[numelements + 1];
                        JRadioButton[] buttons = new JRadioButton[numelements];
                        for (short i = 0; i < numelements; i++) {
                            buttons[i] = new javax.swing.JRadioButton("Risk Group " + (i + currentHighest));
                            group.add(buttons[i]);
                            array[i + 1] = buttons[i];
                        }
                        array[0] = msgString;

                        /* // sbhatia commented this out
                        int sel = -1;
                        while (sel != 0) {
                            sel = JOptionPane.showConfirmDialog(null, array, "", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE);
                        }*/
                        /*buttons[0].setSelected(true); // sbhatia added this line
                        scanButtons:
                        for (short i = 1; i < numelements; i++) {
                            if (buttons[i].isSelected()) {
                                relayRiskGroup((short) (currentHighest + 1));
                                break scanButtons;
                            }
                        }
                        currentHighest = _partDatum._riskGroup;

                    } else {
                        firsthigher = true;
                    }
                }
            }

            //If it's a basic Part, the risk group is whatever the algorithm said
        } else {
            relayRiskGroup(rg);
        }

        //Check it's features to see if any are higher RG
        for (Annotation an : this.getSeq().getAnnotations()) {
            Feature afeat = an.getFeature();
            if (afeat == null) {
                continue;
            }
            relayRiskGroup(afeat.getRiskGroup());
        }
    }

    private void relayRiskGroup(short value) {
        if (value > _partDatum._riskGroup) {
            _partDatum._riskGroup = value;
            System.out.println("Setting risk group to " + _partDatum._riskGroup);
            setChanged(RefreshEvent.Condition.RISK_GROUP_CHANGED);
        } else {
            fireData(new RefreshEvent(this, RefreshEvent.Condition.RISK_GROUP_CHANGED));
        }
    }*/

    public static Part retrieveByName(String name) {
        //query connection for one part whose name contains the provided string    
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public static Part retrieveByExactName(String name) {
        //query connection for one part whose name is the provided string
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @AssertTrue
    private boolean checkDBConstraints(){
        //name must be unique
        //format(by name) + sequence(by content)combination should be unique
        return true;
    }

    public static enum PartType {
        BASIC, COMPOSITE
    };
    
    public static enum PartFunction {
        //XXX: composite is not a function
        CDS, RBS, PROMOTER, TERMINATOR, COMPOSITE;
    }
}
