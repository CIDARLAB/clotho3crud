package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
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

    public Parameter createParameter(double value, Variable variable) {
        Parameter parameter = new Parameter(value, variable);
        addParameter(parameter);
        return parameter;
    }
    
    public void addParameter(Parameter parameter) {
    	if (parameters == null) {
    		parameters = new HashSet<Parameter>();
    	}
    	parameters.add(parameter);
    }
    
    public void addPart(Part part) {
    	if (parts == null) {
    		parts = new HashSet<Part>();
    	}
    	parts.add(part);
    }
    
    public void addPolynucleotide(Polynucleotide polynucleotide) {
    	if (polynucleotides == null) {
    		polynucleotides = new HashSet<Polynucleotide>();
    	}
    	polynucleotides.add(polynucleotide);
    }
    
    public void addStrain(Strain strain) {
    	if (strains == null) {
    		strains = new HashSet<Strain>();
    	}
    	strains.add(strain);
    }
    
    public void addMedium(Medium medium) {
    	if (media == null) {
    		media = new HashSet<Medium>();
    	}
    	media.add(medium);
    }
    
    public void addSubDesign(BioDesign subDesign) {
    	if (subDesigns == null) {
    		subDesigns = new HashSet<BioDesign>();
    	}
    	subDesigns.add(subDesign);
    }

}
