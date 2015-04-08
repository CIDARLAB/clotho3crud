package org.clothocad.model;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import lombok.ToString;
import lombok.extern.slf4j.Slf4j;

/**
 *
 * @author jcanderson
 */
@ToString(callSuper=true, includeFieldNames=true)
@NoArgsConstructor
@Slf4j
public class Oligo extends ObjBase {

    @Setter
    @Getter
    @Reference
    private Person author;

    @Getter
    @Setter
    private String sequence;

    @Getter
    @Setter
    private String description;

    protected Oligo(String name, String desc, String sequence, Person author) {
        super(name); 
        this.description = desc;
        this.author = author;
        this.sequence = sequence;
    }

}
