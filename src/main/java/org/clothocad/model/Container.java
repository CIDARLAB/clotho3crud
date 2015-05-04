package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Container extends SharableObjBase {

    @Getter
    @Setter
    @Reference
    protected Sample sample;

    @Getter
    @Setter
    @Reference
    protected Plate plate;

    public Container(String name, Person author) {
        super(name, author);
    }

    public Container(String name, String description, Person author) {
        super(name, author, description);
    }

}
