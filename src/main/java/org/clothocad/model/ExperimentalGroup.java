package org.clothocad.model;

import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

@NoArgsConstructor
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
    
    public void addSample(Sample sample) {
    	if (samples == null) {
    		samples = new HashSet<Sample>();
    	}
    	samples.add(sample);
    }
    
    public void addSampleData(SampleData sampleDatum) {
    	if (sampleData == null) {
    		sampleData = new HashSet<SampleData>();
    	}
    	sampleData.add(sampleDatum);
    }

}
