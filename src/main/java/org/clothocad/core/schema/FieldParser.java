/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.persistence.annotations.Reference;
import org.objectweb.asm.AnnotationVisitor;
import org.objectweb.asm.FieldVisitor;
import static org.objectweb.asm.Opcodes.*;
import org.objectweb.asm.Type;
/**
 *
 * @author spaige
 */
public class FieldParser extends FieldVisitor{

    private ClothoField target;
    public FieldParser(ClothoField target) {
        super(ASM4);
        this.target = target;
    }

    @Override
    public AnnotationVisitor visitAnnotation(String desc, boolean visible) {
        if (desc.equals(Type.getType(Reference.class))) target.setReference(true);
        return null;
    }
    
    
    
}
