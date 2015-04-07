package org.clothocad.model;

import java.util.HashMap;
import java.util.List;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

/**
*
* @author Nicholas Roehner
*/
public class Plate extends SharableObjBase {
	
	@Getter
	@Setter
	@ReferenceCollection
	protected List<Container> containers;
	
	@Getter
	@Setter
	protected HashMap<Container, List<Integer>> containerCoordinates;
	
	public Plate(String name, Person author) {
		super(name, author);
	}
	
	public Plate(String name, String description, Person author) {
		super(name, author, description);
	}

}
