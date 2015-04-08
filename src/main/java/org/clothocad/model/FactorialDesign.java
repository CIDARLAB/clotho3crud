package org.clothocad.model;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

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
public class FactorialDesign extends ExperimentalDesign {

    @NotNull
    @Size(min=1)
    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Factor> factors = new HashSet<Factor>();

    public FactorialDesign(String name, Set<Variable> responseVariables, Set<Variable> controlledVariables,
            Person author) {
        super(name, responseVariables, controlledVariables, author);
    }

    public FactorialDesign(String name, String description, Set<Variable> responseVariables,
            Set<Variable> controlledVariables, Person author) {
        super(name, description, responseVariables, controlledVariables, author);
    }

    public Factor createFactor(String name, Variable variable) {
        if (factors == null) {
            factors = new HashSet<Factor>();
        }
        Factor factor = new Factor(name, variable);
        factors.add(factor);
        return factor;
    }

}
