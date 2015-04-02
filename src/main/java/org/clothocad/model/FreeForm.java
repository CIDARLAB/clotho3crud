/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.List;
import org.clothocad.model.Format;
import org.clothocad.model.NucSeq;
import org.clothocad.model.Part;

/**
 *
 * @author spaige
 */
public class FreeForm implements Format {

	@Override
    public boolean checkPart(Part p) {
        return true;
    }

	@Override
    public boolean checkComposite(List<Part> composition, Object additionalRequirements) {
        return true;
    }

	@Override
    public NucSeq generateCompositeSequence(List<Part> composition, Object additionalRequirements) {
        //XXX: dummy implementation
        StringBuilder builder = new StringBuilder();
        for (Part part : composition){
            builder.append(part.getSequence().getSeq());
        }
        return new NucSeq(builder.toString());
    }
    
}
