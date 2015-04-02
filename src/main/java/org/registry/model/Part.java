package org.registry.model;



import java.util.ArrayList;
import java.util.List;

import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

@Data()
@NoArgsConstructor
public class Part extends SharableObjBase {

	private List<Feature> features = new ArrayList<Feature>();
	private List<Subscar> specified_subscars = new ArrayList<Subscar>();
	private List<Subpart> specified_subparts = new ArrayList<Subpart>();
	private List<Subpart> deep_subparts = new ArrayList<Subpart>();
	private List<String> sequences = new ArrayList<String>();
	private List<String> twins = new ArrayList<String>();
	

	private String release_status, references, part_nickname, parameters,
		part_url, part_type, sample_status, part_results, samples,
		part_short_name, part_rating, part_id, part_short_desc, categories, groups,
		part_name, part_author, part_entered;
}