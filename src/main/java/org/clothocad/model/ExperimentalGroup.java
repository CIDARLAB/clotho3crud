package org.clothocad.model;

import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.Setter;

public class ExperimentalGroup {
	
	@NotNull
	@Size(min=1)
	@Getter
	@Setter
	@ReferenceCollection
	protected Set<Sample> samples;
	
	@Getter
	@Setter
	@ReferenceCollection
	protected Set<SampleData> sampleData;
	
	@Getter
	@Setter
	@Reference
	protected ExperimentalCondition experimentalCondition;
	
	protected ExperimentalGroup(Set<Sample> samples) {
		this.samples = samples;
	}
	
	

}
