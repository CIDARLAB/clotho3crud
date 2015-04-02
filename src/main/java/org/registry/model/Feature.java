package org.registry.model;



import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

@Data()
@NoArgsConstructor
public class Feature extends SharableObjBase {

	private boolean direction;
	private String title, type;
	private int startpos, endpos;

}