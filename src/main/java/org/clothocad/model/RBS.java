package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;

import lombok.Data;
import lombok.NoArgsConstructor;

@Data()
@NoArgsConstructor
public class RBS extends SharableObjBase {

    // native_context --> nativeContext
    // native_host --> nativeHost
    private String utr5, nativeContext, nativeHost, sequence, weblink, isIntergenic;

}
