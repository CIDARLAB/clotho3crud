/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.List;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.persistence.annotations.Add;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

/**
 *
 * @author spaige
 */
@NoArgsConstructor
@Add(name="sequence", provider="getSequence")
public class CompositePart extends Part {

    @Getter
    @Setter
    @ReferenceCollection
    private List<Part> composition;
        
    private static final PartType partType= PartType.COMPOSITE;

    public CompositePart(List<Part> composition, Object additionalRequirements, Format f, Person author, String name, String shortdescription) {
        super(name, shortdescription, f, author);
        if (!f.checkComposite(composition, additionalRequirements)) {
            System.out.println("generateComposite: Doesn't obey format, return null");
        }
        setComposition(composition);
        setType(PartFunction.COMPOSITE);
    }

    @Override
    public boolean checkFormat() {
          return getFormat().checkComposite(this.composition, null);
    }

    @Override
    public PartType getPartType() {
        return partType;
    }

    @Override
    public NucSeq getSequence() {
                    return getFormat().generateCompositeSequence(composition, null);
    }
}
