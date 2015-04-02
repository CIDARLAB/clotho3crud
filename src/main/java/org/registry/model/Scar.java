package org.registry.model;



import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

@Data()
@NoArgsConstructor
public class Scar extends SharableObjBase implements Subscar {

	private String scar_sequence, scar_nickname, scar_comments, scar_standard,
		scar_type, name;

}