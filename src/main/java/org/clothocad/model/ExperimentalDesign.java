package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.ReferenceCollection;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

/**
*
* @author Nicholas Roehner
*/
public class ExperimentalDesign extends SharableObjBase {

    @NotNull
    @Size(min=1)
    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Variable> responseVariables;

    @NotNull
    @Size(min=1)
    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Variable> controlledVariables;

    @Getter
    protected Set<ExperimentalCondition> experimentalConditions;

    @Getter
    @Setter
    @Reference
    protected BioDesign bioDesign;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<ExperimentalDesign> subDesigns;

    @Getter
    @Setter
    @Reference
    protected ExperimentalDesign parentDesign;

    public ExperimentalDesign(String name, Set<Variable> responseVariables, Set<Variable> controlledVariables,
            Person author) {
        super(name, author);
        this.responseVariables = responseVariables;
        this.controlledVariables = controlledVariables;
    }

    public ExperimentalDesign(String name, String description, Set<Variable> responseVariables,
            Set<Variable> controlledVariables, Person author) {
        super(name, author, description);
        this.responseVariables = responseVariables;
        this.controlledVariables = controlledVariables;
    }

    public ExperimentalCondition createExperimentalCondition() {
        if (experimentalConditions == null) {
            experimentalConditions = new HashSet<ExperimentalCondition>();
        }
        ExperimentalCondition experimentalCondition = new ExperimentalCondition();
        experimentalConditions.add(experimentalCondition);
        return experimentalCondition;
    }

}
