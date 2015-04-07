package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;

public abstract class Derivation extends SharableObjBase {
	
	public Derivation(String name, Person author) {
		super(name, author);
	}
	
	public Derivation(String name, String description, Person author) {
		super(name, author);
	}
	
	public abstract Parameter deriveParameter(SampleData sampleData);
	
}
