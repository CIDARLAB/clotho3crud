package org.clothocad.model;

import lombok.NoArgsConstructor;

import org.clothocad.core.datums.SharableObjBase;

@NoArgsConstructor
public class Instrument extends SharableObjBase {

    public Instrument(String name, Person author) {
        super(name, author);
    }

    public Instrument(String name, String description, Person author) {
        super(name, author, description);
    }

}
