package org.clothocad.model;

import java.util.ArrayList;
import java.util.List;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

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
	
	protected Assembly() {
		
	}
	
	public Assembly createSubAssembly() {
		if (subAssemblies == null) {
			subAssemblies = new ArrayList<Assembly>();
    	}
		Assembly subAssembly = new Assembly();
		subAssemblies.add(subAssembly);
		return subAssembly;
	}
}
