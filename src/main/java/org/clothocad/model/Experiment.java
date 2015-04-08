package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.ReferenceCollection;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;

/**
*
* @author Nicholas Roehner
*/
public class Experiment extends SharableObjBase {

    @NotNull
    @Getter
    @Setter
    @Reference
    protected ExperimentalDesign experimentalDesign;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Experiment> subExperiments;

    @Getter
    @Setter
    @Reference
    protected Experiment parentExperiment;

    @Getter
    protected Set<ExperimentalGroup> experimentalGroups;

    public Experiment(String name, ExperimentalDesign experimentalDesign, Person author) {
        super(name, author);
        this.experimentalDesign = experimentalDesign;
    }

    public Experiment(String name, String description, ExperimentalDesign experimentalDesign, Person author) {
        super(name, author, description);
        this.experimentalDesign = experimentalDesign;
    }

    public ExperimentalGroup createExperimentalGroup(Set<Sample> samples) {
        if (experimentalGroups == null) {
            experimentalGroups = new HashSet<ExperimentalGroup>();
        }
        ExperimentalGroup experimentalGroup = new ExperimentalGroup(samples);
        experimentalGroups.add(experimentalGroup);
        return experimentalGroup;
    }

}
