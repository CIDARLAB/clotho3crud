package org.clothocad.model;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.ArrayList;
import java.util.List;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Assembly {

    @Getter
    @Setter
    @ReferenceCollection
    protected List<Part> parts;

    @Getter
    protected List<Assembly> subAssemblies;
    
    public Assembly createSubAssembly() {
        Assembly subAssembly = new Assembly();
       	addSubAssembly(subAssembly);
        return subAssembly;
    }
    
    public void addPart(Part part) {
    	if (parts == null) {
    		parts = new ArrayList<Part>();
    	}
    	parts.add(part);
    }
    
    public void addSubAssembly(Assembly subAssembly) {
    	if (subAssemblies == null) {
    		subAssemblies = new ArrayList<Assembly>();
    	}
    	subAssemblies.add(subAssembly);
    }
}
