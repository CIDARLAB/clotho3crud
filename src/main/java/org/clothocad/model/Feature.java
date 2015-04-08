package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import javax.validation.constraints.NotNull;

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
