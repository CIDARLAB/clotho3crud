package org.clothocad.model;

import lombok.Getter;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

/**
*
* @author Nicholas Roehner
*/
public class ExperimentalCondition {

    @NotNull
    @Size(min=1)
    @Getter
    protected Set<Parameter> parameters;

    protected ExperimentalCondition() {}

    public Parameter createParameter(double value, Variable variable, Units units) {
        if (parameters == null) {
            parameters = new HashSet<Parameter>();
        }
        Parameter parameter = new Parameter(value, variable, units);
        parameters.add(parameter);
        return parameter;
    }

}
