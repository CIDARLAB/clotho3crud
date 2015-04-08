package org.clothocad.model;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.Setter;

import java.util.ArrayList;
import java.util.List;

/**
*
* @author Nicholas Roehner
*/
public class Assembly {

    @Getter
    @Setter
    @ReferenceCollection
    protected List<Part> parts;

    @Getter
    protected List<Assembly> subAssemblies;

    protected Assembly() {}
    
    public Assembly createSubAssembly() {
        if (subAssemblies == null) {
            subAssemblies = new ArrayList<Assembly>();
        }
        Assembly subAssembly = new Assembly();
        subAssemblies.add(subAssembly);
        return subAssembly;
    }
}
