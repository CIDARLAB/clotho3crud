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

import javax.validation.constraints.NotNull;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;

/**
 *
 * @author J. Christopher Anderson
 * @author Nicholas Roehner
 */
@NoArgsConstructor
public class Feature extends SharableObjBase {
	
	@Setter
    @Getter
    @Reference
    protected Sequence sequence;
    
    @Setter
    @Getter
    protected String genbankId, swissProtId, PDBId;
    
    @Getter
    protected short riskGroup;
   
    @NotNull
    @Setter
    @Getter
    protected FeatureRole role;
    
    @Getter
	@Setter
	@Reference
	protected Feature parentFeature;
    
    /**
     * Constructor of a new Feature
     * @param name
     * @param role
     * @param author
     */
    public Feature(String name, FeatureRole role, Person author) {
        super(name, author);
        this.role = role;
    }
    
    /**
     * Constructor of a new Feature
     * @param name
     * @param description
     * @param role
     * @param author
     */
    public Feature(String name, String description, FeatureRole role, Person author) {
        super(name, author, description);
        this.role = role;
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
    
    // Feel free to add more of these
    public static enum FeatureRole {
    	PROMOTER, CDS, RBS, TERMINATOR, SCAR, SPACER, RIBOZYME;
    }

}
