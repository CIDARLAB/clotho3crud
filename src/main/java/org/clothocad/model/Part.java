package org.clothocad.model;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.validation.constraints.NotNull;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.model.Feature.FeatureRole;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
 * @author Nicholas Roehner
 */
@NoArgsConstructor
public class Part extends SharableObjBase {
	
	@Getter
	@Reference
	protected Format format;
	
	@Getter
	protected List<Assembly> assemblies;
	
	@NotNull
	@Getter
	@Setter
	@Reference
	protected Sequence sequence;
	
	@Getter
    @Setter
    protected boolean isForwardOrientation;
	
	@Getter
	@Setter
	@Reference
	protected Part parentPart;
    
    public Part(String name, String description, Sequence sequence, Person author){
        super(name, author, description); 
        this.sequence = sequence;
    }
    
    public Part(String name, Sequence sequence, Person author){
        super(name, author); 
        this.sequence = sequence;
    }

    /**
     * Change the Format of the Part
     * @param format new Format for the Part
     */
    public void setFormat(Format format) {
        if (format.checkPart(this)) {
        	this.format = format;
        }
    }
    
    public List<FeatureRole> getRoles() {
    	List<FeatureRole> roles = new LinkedList<FeatureRole>();
    	for (Annotation annotation : sequence.getAnnotations()) {
    		Feature feature = annotation.getFeature();
    		if (feature != null) {
    			roles.add(feature.getRole());
    		}
    	}
    	return roles;
    }
    
    public Assembly createAssembly() {
    	if (assemblies == null) {
    		assemblies = new ArrayList<Assembly>();
    	}
    	Assembly assembly = new Assembly();
    	assemblies.add(assembly);
    	return assembly;
    }
    
}
