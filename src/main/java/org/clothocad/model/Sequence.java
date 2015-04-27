package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.ArrayList;
import java.util.List;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Pattern;

/**
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Sequence extends SharableObjBase {

    @NotNull
    @Getter
    @Setter
    @Pattern(regexp="[ATUCGRYKMSWBDHVN]*", flags={Pattern.Flag.CASE_INSENSITIVE})
    protected String sequence;

    @Getter
    @Setter
    protected List<Annotation> annotations;

    @Getter
    @Setter
    @Reference
    protected Sequence parentSequence;

    public Sequence(String name, String sequence, Person author) {
        super(name, author);
        this.sequence = sequence;
    }

    public Sequence(String name, String description, String sequence, Person author) {
        super(name, author, description);
        this.sequence = sequence;
    }

    public Annotation createAnnotation(String name, int start, int end, boolean isForwardStrand,
            Person author) {
        if (annotations == null) {
            annotations = new ArrayList<Annotation>();
        }
        Annotation annotation = new Annotation(name, start, end, isForwardStrand, author);
        annotations.add(annotation);
        return annotation;
    }

    public Annotation createAnnotation(String name, String description, int start, int end,
            boolean isForwardStrand, Person author) {
        if (annotations == null) {
            annotations = new ArrayList<Annotation>();
        }
        Annotation annotation = new Annotation(name, description, start, end, isForwardStrand, author);
        annotations.add(annotation);
        return annotation;
    }

}
