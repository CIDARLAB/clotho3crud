package org.clothocad.model;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Factor extends ObjBase {

    @NotNull
    @Size(min=2)
    @Getter
    @ReferenceCollection
    protected Set<Level> levels;

    @NotNull
    @Getter
    @Setter
    @Reference
    protected Variable variable;

    @Getter
    @Setter
    @Reference
    protected BioDesign bioDesign;

    @Getter
    @Setter
    protected String description;

    protected Factor(String name, Variable variable) {
        super(name);
        this.variable = variable;
    }

    public Level createLevel(String name) {
        Level level = new Level(name);
        addLevel(level);
        return level;
    }
    
    public void addLevel(Level level) {
    	if (levels == null) {
    		levels = new HashSet<Level>();
    	}
    	levels.add(level);
    }

}
