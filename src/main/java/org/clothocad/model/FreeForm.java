package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;

import lombok.NoArgsConstructor;

import java.util.List;

/**
 *
 * @author spaige
 * @author Nicholas Roehner
 */
@NoArgsConstructor
public class FreeForm extends SharableObjBase implements Format {

    public FreeForm(Person author) {
        super("FreeForm", author);
    }

    public FreeForm(String description, Person author) {
        super("FreeForm", author, description);
    }

    public boolean checkPart(Part p) {
        return true;
    }

    public boolean checkComposite(List<Part> subParts) {
        return true;
    }

    public Part generateCompositePart(String name, List<Part> subParts, Person author) {
        Sequence compositeSeq = buildCompositeSequence(subParts, author);
        Part compositePart = new Part(name, compositeSeq, author);
        compositePart.setFormat(this);
        assembleCompositePart(compositePart, subParts);
        return compositePart;
    }

    public Part generateCompositePart(String name, String description, List<Part> subParts, Person author) {
        Sequence compositeSeq = buildCompositeSequence(subParts, author);
        Part compositePart = new Part(name, description, compositeSeq, author);
        compositePart.setFormat(this);
        assembleCompositePart(compositePart, subParts);
        return compositePart;
    }

    private Sequence buildCompositeSequence(List<Part> subParts, Person author) {
        StringBuilder builder = new StringBuilder();
        for (Part part : subParts) {
            builder.append(part.getSequence().getSequence());
        }
        return new Sequence("seq", builder.toString(), author);
    }

    private void assembleCompositePart(Part compositePart, List<Part> subParts) {
        Assembly assembly = compositePart.createAssembly();
        assembly.setParts(subParts);
    }

}
