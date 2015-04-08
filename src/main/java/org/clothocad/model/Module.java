package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.Setter;

import java.util.Set;

import javax.validation.constraints.NotNull;

/**
*
* @author Nicholas Roehner
*/
public abstract class Module extends SharableObjBase {

    @NotNull
    @Getter
    @Setter
    protected ModuleRole role;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Influence> influences;

    @Getter
    @Setter
    @Reference
    protected Module parentModule;

    public Module(String name, ModuleRole role, Person author) {
        super(name, author);
        this.role = role;
    }

    public Module(String name, String description, ModuleRole role, Person author) {
        super(name, author, description);
        this.role = role;
    }

    // Feel free to add more of these
    public static enum ModuleRole {
        TRANSCRIPTION, TRANSLATION, EXPRESSION;
    }

}
