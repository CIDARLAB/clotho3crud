package org.clothocad.model;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

public abstract class SampleData extends SharableObjBase {
	
	@NotNull
	@Getter
	@Setter
	@Reference
	protected Sample sample;
	
	@NotNull
	@Getter
	@Setter
	@Reference
	protected Instrument instrument;
	
	@NotNull
	@Size(min=1)
	@Getter
	protected Set<Parameter> parameters;
	
	public SampleData(String name, Sample sample, Instrument instrument, Set<Variable> responseVariables, 
			Person author) {
		super(name, author);
		this.sample = sample;
		this.instrument = instrument;
	}
	
	public SampleData(String name, String description, Sample sample, Instrument instrument,
			Set<Variable> responseVariables, Person author) {
		super(name, author, description);
		this.sample = sample;
		this.instrument = instrument;
	}
	
	public Parameter createParameter(double value, Variable variable, Units units) {
		if (parameters == null) {
			parameters = new HashSet<Parameter>();
		}
		Parameter parameter = new Parameter(value, variable, units);
		parameters.add(parameter);
		return parameter;
	}
}
