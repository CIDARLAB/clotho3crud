package org.clothocad.model;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.Setter;

import java.util.List;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

/**
*
* @author Nicholas Roehner
*/
public class CompositeModule extends Module {

    @NotNull
    @Size(min=2)
    @Getter
    @Setter
    @ReferenceCollection
    protected List<Module> subModules;

    public CompositeModule(String name, ModuleRole role, List<Module> subModules, Person author) {
        super(name, role, author);
        this.subModules = subModules;
    }

    public CompositeModule(String name, String description, ModuleRole role, List<Module> subModules,
            Person author) {
        super(name, description, role, author);
        this.subModules = subModules;
    }

}
