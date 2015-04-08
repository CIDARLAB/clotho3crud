package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

/**
*
* @author Nicholas Roehner
*/
public class BioDesign extends SharableObjBase {

    @Getter
    @Setter
    @Reference
    protected Module module;

    @Getter
    protected Set<Parameter> parameters;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Part> parts;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Polynucleotide> polynucleotides;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Strain> strains;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Medium> media;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<BioDesign> subDesigns;

    @Getter
    @Setter
    @Reference
    protected BioDesign parentDesign;

    public BioDesign(String name, String description, Person author) {
        super(name, author, description);
    }

    public BioDesign(String name, Person author) {
        super(name, author);
    }

    public Parameter createParameter(double value, Variable variable, Units units) {
        if (parameters == null) {
            parameters = new HashSet<Parameter>();
        }
        Parameter parameter = new Parameter(value, variable, units);
        parameters.add(parameter);
        return parameter;
    }

}
