package org.clothocad.model;

import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.Setter;

import javax.validation.constraints.NotNull;

/**
*
* @author Nicholas Roehner
*/
public class Parameter {

    @NotNull
    @Getter
    @Setter
    protected double value;

    @NotNull
    @Getter
    @Setter
    @Reference
    protected Variable variable;

    @NotNull
    @Getter
    @Setter
    @Reference
    protected Units units;

    @Getter
    @Setter
    @Reference
    protected Derivation derivation;

    protected Parameter(double value, Variable variable, Units units) {
        this.value = value;
        this.variable = variable;
        this.units = units;
    }

}
