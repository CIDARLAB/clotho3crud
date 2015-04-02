package org.clothocad.model;

import java.util.ArrayList;
import java.util.List;

import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

@Data()
@NoArgsConstructor
public class SinglePathway extends SharableObjBase {

	private List<Polypeptide> intermediates = new ArrayList<Polypeptide>();
	private List<Reaction> reactions = new ArrayList<Reaction>();

	private String version, target;

}