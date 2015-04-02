package org.clothocad.model;

import javax.validation.constraints.Pattern;
import lombok.Getter;
import lombok.Setter;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
public class SimpleSequence extends ObjBase {
    
    public SimpleSequence(String name, String sequence){
        super(name);
        this.sequence = sequence;
    }
    @Getter
    @Setter
    @Pattern(regexp="[ATUCGRYKMSWBDHVN]*", flags={Pattern.Flag.CASE_INSENSITIVE})
    String sequence;
    
    
    
}
