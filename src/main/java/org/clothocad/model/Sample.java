package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import javax.validation.constraints.NotNull;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Sample extends SharableObjBase {

    @Getter
    @Setter
    @Reference
    protected BioDesign bioDesign;

    @NotNull
    @Getter
    @Setter
    @Reference
    protected Container container;

    @Getter
    @Setter
    @Reference
    protected Sample parentSample;

    public Sample(String name, Container container, Person author) {
        super(name, author);
        this.container = container;
    }

    public Sample(String name, String description, Container container, Person author) {
        super(name, author, description);
        this.container = container;
    }

}
