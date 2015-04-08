package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;

import lombok.Data;
import lombok.NoArgsConstructor;

import java.util.ArrayList;
import java.util.List;

@Data()
@NoArgsConstructor
public class SinglePathway extends SharableObjBase {

    private List<Polypeptide> intermediates = new ArrayList<Polypeptide>();
    private List<Reaction> reactions = new ArrayList<Reaction>();

    private String version, target;

}
