package org.clothocad.model;

import javax.validation.constraints.NotNull;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

/**
*
* @author Nicholas Roehner
*/
public class Influence extends SharableObjBase {
	
	@NotNull
	@Getter
	@Setter
	@Reference
	protected Feature influencingFeature;
	
	@NotNull
	@Getter
	@Setter
	@Reference
	protected Feature influencedFeature;
	
	@NotNull
	@Getter
	@Setter
	protected InfluenceType type;
	
	@Getter
	@Setter
	@Reference
	protected Influence parentInfluence;

	public Influence(String name, Feature influencingFeature, Feature influencedFeature, 
			InfluenceType type, Person author) {
		super(name, author);
		this.influencingFeature = influencingFeature;
		this.influencedFeature = influencedFeature;
		this.type = type;
	}
	
	public Influence(String name, String description, Feature influencingFeature, Feature influencedFeature, 
			InfluenceType type, Person author) {
		super(name, author, description);
		this.influencingFeature = influencingFeature;
		this.influencedFeature = influencedFeature;
		this.type = type;
	}
	
	// Feel free to add more of these
    public static enum InfluenceType {
    	REPRESSION, ACTIVATION;
    }
}
