package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.ReferenceCollection;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
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
        ExperimentalCondition experimentalCondition = new ExperimentalCondition();
        addExperimentalCondition(experimentalCondition);
        return experimentalCondition;
    }
    
    public void addExperimentalCondition(ExperimentalCondition experimentalCondition) {
    	if (experimentalConditions == null) {
    		experimentalConditions = new HashSet<ExperimentalCondition>();
    	}
    	experimentalConditions.add(experimentalCondition);
    }
    
    public void addSubDesign(ExperimentalDesign subDesign) {
    	if (subDesigns == null) {
    		subDesigns = new HashSet<ExperimentalDesign>();
    	}
    	subDesigns.add(subDesign);
    }
    
    public void addResponseVariable(Variable responseVariable) {
    	if (responseVariables == null) {
    		responseVariables = new HashSet<Variable>();
    	}
    	responseVariables.add(responseVariable);
    }
    
    public void addControlledVariable(Variable controlledVariable) {
    	if (controlledVariables == null) {
    		controlledVariables = new HashSet<Variable>();
    	}
    	controlledVariables.add(controlledVariable);
    }

}
