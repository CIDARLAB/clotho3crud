/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import java.util.HashSet;
import org.clothocad.core.datums.Argument;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.objectweb.asm.ClassReader;
import org.objectweb.asm.ClassVisitor;
import org.objectweb.asm.FieldVisitor;
import org.objectweb.asm.MethodVisitor;
import static org.objectweb.asm.Opcodes.*;
import org.objectweb.asm.Type;

/**
 *
 * @author spaige
 */
public class ClassParser extends ClassVisitor{
    private Schema target;
    
    public ClassParser(Schema target){
        super(ASM4);
        target.methods = new HashSet<>();
        target.fields = new HashSet<>();
    }

    @Override
    public MethodVisitor visitMethod(int access, String name, String desc, String signature, String[] exceptions) {
        String description = "";
        
        Class[] input = typesToClasses(Type.getArgumentTypes(desc));
        Argument[] arguments = new Argument[input.length];
        for (int i = 0; i<input.length; i ++){
            //TODO: real argument names
            arguments[i] = new Argument("", input[i]);
        }
        
        Class output = Type.getReturnType(desc).getClass();
        
        
        Function method = new Function(name, arguments, output, "", target.getLanguage());
        target.methods.add(method);
        return null;
    }

    @Override
    public FieldVisitor visitField(int access, String name, String desc, String signature, Object value) {
        String example = "";
        String description = "";
        Access acc = opcodeToAccess(access);
        Boolean reference = false;
        Class type = Type.getType(desc).getClass();
        
        ClothoField field = new ClothoField(name, type, example, description, reference, acc);
        
        target.fields.add(field);
        return new FieldParser(field);
    }
    
    public static Access opcodeToAccess(int opcode){
        if ((opcode & ACC_PUBLIC) != 0) return Access.PUBLIC;
        return Access.PRIVATE;
    }
    
    public static Class[] typesToClasses(Type[] types){
        Class[] result = new Class[types.length];
        for (int i = 0; i < types.length; i++){
            result[i] = types.getClass();
        }
        return result;
    }
    
    
}
