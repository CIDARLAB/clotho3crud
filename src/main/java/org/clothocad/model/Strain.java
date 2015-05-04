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
public class Strain extends SharableObjBase {

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Polynucleotide> polynucleotides;

    @Getter
    @Setter
    @Reference
    protected Strain parentStrain;

    public Strain(String name, Person author) {
        super(name, author);
    }

    public Strain(String name, String description, Person author) {
        super(name, author, description);
    }
    
    public void addPolynucleotide(Polynucleotide polynucleotide) {
    	if (polynucleotides == null) {
    		polynucleotides = new HashSet<Polynucleotide>();
    	}
    	polynucleotides.add(polynucleotide);
    }

}
