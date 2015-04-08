package org.clothocad.model;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

import lombok.Data;
import lombok.NoArgsConstructor;

import java.util.ArrayList;
import java.util.List;

@Data()
@NoArgsConstructor
public class Enzyme extends SharableObjBase {

    private List<Polypeptide> polypeptides = new ArrayList<Polypeptide>();

    private String version;

}
