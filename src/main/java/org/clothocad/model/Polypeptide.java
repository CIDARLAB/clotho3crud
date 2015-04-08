package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;

import lombok.Data;
import lombok.NoArgsConstructor;

@Data()
@NoArgsConstructor
public class Polypeptide extends SharableObjBase {

    private String version, source, sequence;

}
