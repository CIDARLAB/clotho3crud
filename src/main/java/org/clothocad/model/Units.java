package org.clothocad.model;

import lombok.NoArgsConstructor;

import org.clothocad.core.datums.SharableObjBase;

@NoArgsConstructor
public class Units extends SharableObjBase {

    public Units(String name, Person author) {
        super(name, author);
    }

    public Units(String name, String description, Person author) {
        super(name, author, description);
    }

}
