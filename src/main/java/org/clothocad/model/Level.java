package org.clothocad.model;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.Setter;

import javax.validation.constraints.NotNull;

/**
*
* @author Nicholas Roehner
*/
public class Level extends ObjBase {

    @NotNull
    @Getter
    protected Parameter parameter;

    @Getter
    @Setter
    @Reference
    protected BioDesign bioDesign;

    @Getter
    @Setter
    protected String description;

    protected Level(String name) {
        super(name);
    }

    public Parameter createParameter(double value, Variable variable, Units units) {
        parameter = new Parameter(value, variable, units);
        return parameter;
    }

}
