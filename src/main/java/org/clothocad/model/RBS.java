package org.clothocad.model;

import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

@Data()
@NoArgsConstructor
public class RBS extends SharableObjBase {

	// native_context --> nativeContext
	// native_host --> nativeHost
	private String utr5, nativeContext, nativeHost, sequence, weblink, isIntergenic;

}