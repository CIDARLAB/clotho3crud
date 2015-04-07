package org.clothocad.model;

import java.util.ArrayList;
import java.util.List;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

/**
*
* @author Nicholas Roehner
*/
public class Factor extends ObjBase {
	
	@NotNull
	@Size(min=2)
	@Getter
	@ReferenceCollection
	protected List<Level> levels;
	
	@NotNull
	@Getter
	@Setter
	@Reference
	protected Variable variable;
	
	@Getter
	@Setter
	@Reference
	protected BioDesign bioDesign;
	
	@Getter
	@Setter
	protected String description;
	
	protected Factor(String name, Variable variable) {
		super(name);
		this.variable = variable;
	}
	
	public Level createLevel(String name) {
		if (levels == null) {
			levels = new ArrayList<Level>();
		}
		Level level = new Level(name);
		levels.add(level);
		return level;
	}

}
