package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;

import lombok.Data;
import lombok.NoArgsConstructor;

import java.util.ArrayList;
import java.util.List;

@Data()
@NoArgsConstructor
public class Reaction extends SharableObjBase {

    private List<Enzyme> enzymes = new ArrayList<Enzyme>();
    private List<Polypeptide> products = new ArrayList<Polypeptide>();

    private String version;

}
