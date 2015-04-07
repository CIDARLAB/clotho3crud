package org.clothocad.model;

import java.util.List;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.Setter;

/**
*
* @author Nicholas Roehner
*/
public class BasicModule extends Module {
	
	@NotNull
	@Size(min=1)
	@Getter
	@Setter
	@ReferenceCollection
	protected List<Feature> features;
	
	public BasicModule(String name, ModuleRole role, List<Feature> features, Person author) {
		super(name, role, author);
		this.features = features;
	}
	
	public BasicModule(String name, String description, ModuleRole role, List<Feature> features, Person author) {
		super(name, description, role, author);
		this.features = features;
	}
	
}
