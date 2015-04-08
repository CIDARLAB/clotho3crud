package org.clothocad.model;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

import java.util.List;

import lombok.Data;
import lombok.NoArgsConstructor;
//import org.clothocad.core.persistence.annotations.Reference;

@Data()
@NoArgsConstructor
public class Highlight extends ObjBase {

    //@Getter
    private int start, end;
    //@Getter
    private String forColor, revColor, inference, description;
    //@Getter
    private boolean plusStrand;
    //@Getter
    private List<String> notes;

    private SharableObjBase refSeq;

    /**
    *public Highlight( int start, int end,
    *   String name, String forColor, String revColor, String description,
    *   String inference, List<String> notes, boolean plusStrand ) {
    *    
    *    super(name);
    *    this.start = start;
    *    this.end = end;
    *    this.forColor = forColor;
    *    this.revColor = revColor;
    *    this.inference = inference;
    *    this.plusStrand = plusStrand;
    *    this.description = description;
    *    this.notes = notes;
    *    this.plusStrand = plusStrand;
    *}
    **/

}