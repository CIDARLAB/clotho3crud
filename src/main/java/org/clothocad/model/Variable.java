package org.clothocad.model;

import lombok.NoArgsConstructor;

import org.clothocad.core.datums.SharableObjBase;

@NoArgsConstructor
public class Variable extends SharableObjBase {

    public Variable(String name, Person author) {
        super(name, author);
    }

    public Variable(String name, String description, Person author) {
        super(name, author, description);
    }

}
