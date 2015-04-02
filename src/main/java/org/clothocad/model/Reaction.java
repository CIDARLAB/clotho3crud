package org.clothocad.model;

import java.util.ArrayList;
import java.util.List;

import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

@Data()
@NoArgsConstructor
public class Reaction extends SharableObjBase {

	private List<Enzyme> enzymes = new ArrayList<Enzyme>();
	private List<Polypeptide> products = new ArrayList<Polypeptide>();

	private String version;

}