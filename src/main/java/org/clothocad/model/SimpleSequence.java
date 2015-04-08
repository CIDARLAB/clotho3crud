package org.clothocad.model;

import lombok.NoArgsConstructor;

/**
 *
 * @author spaige
 * @author Nicholas Roehner
 */
@NoArgsConstructor
public class SimpleSequence extends Sequence {

    public SimpleSequence(String sequence, Person author) {
        super("Simple Seq", sequence, author);
    }

    public SimpleSequence(String name, String sequence, Person author) {
        super(name, sequence, author);
    }

}
