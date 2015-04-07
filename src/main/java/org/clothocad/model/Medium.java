package org.clothocad.model;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

/**
 * 
 * @author Nicholas Roehner
 */ 
public class Medium extends SharableObjBase {
	
	@Getter
	@Setter
	@Reference
	protected Strain parentMedium;
	
	public Medium(String name, Person author) {
		super(name, author);
	}
	
	public Medium(String name, String description, Person author) {
		super(name, author, description);
	}

}
