/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.model.Person;
import org.objectweb.asm.AnnotationVisitor;
import org.objectweb.asm.ClassWriter;
import org.objectweb.asm.FieldVisitor;
import org.objectweb.asm.MethodVisitor;
import static org.objectweb.asm.Opcodes.*;
import org.objectweb.asm.Type;
import org.objectweb.asm.util.CheckClassAdapter;

/**
 *
 * @author spaige
 */
@Slf4j
public class ClothoSchema extends Schema {

    public ClothoSchema(String name, String description, Person author, Schema superClass, Set<ClothoField> fields) {
        super(name, description, author);
        this.fields = fields == null? new HashSet<ClothoField>() : fields;
        this.superClass = superClass;
    }

    @JsonCreator
    public ClothoSchema(@JsonProperty("name") String name, @JsonProperty("description") String description, 
        @JsonProperty("author") Person author, @JsonProperty("superClass") Schema superClass, 
        @JsonProperty("fields") Set<ClothoField> fields, @JsonProperty("methods") Set<Function> methods,
        @JsonProperty("dependencies") Map<String,ObjectId> dependencies, @JsonProperty("classData") byte[] classData){
        this(name, description, author, superClass, fields);
        this.methods = methods == null? new HashSet<Function>() : methods;
        this.dependencies = dependencies;
        this.classData = classData;
    }
    
    @Override
    @JsonProperty("language")
    public Language getLanguage() {
        return Language.JSONSCHEMA;
    }

    @Override
    public void setSource(String source) {
        throw new UnsupportedOperationException("Cannot set source on JSON-based schemas"); //TODO: should be equivalent to setJSON?
    }

    @Override
    public byte[] getClassData() {
        if (classData == null) {
            classData = generateClassData();
        }
        return classData;
    }

    protected byte[] generateClassData() {
        log.trace("generating class bytecode for schema {}", this.getId());

        ClassWriter cwriter = new ClassWriter(0);
        //TraceClassVisitor tcv = new TraceClassVisitor(cwriter, new PrintWriter(System.out));
        CheckClassAdapter cw = new CheckClassAdapter(cwriter);

        //ClassWriter cw = new ClassWriter(0);
        String superClassName = this.superClass == null ? "org/clothocad/core/datums/ObjBase" : this.superClass.getInternalName();

        cw.visit(V1_7, ACC_PUBLIC, this.getInternalName(), null, superClassName, new String[]{});
        //store original name
        cw.visitField(ACC_PUBLIC + ACC_FINAL + ACC_STATIC, "SCHEMA_NAME", Type.getType(String.class).getDescriptor(), null, this.getName()).visitEnd();

        //no-args constructor
        log.trace("Creating no-args constructor.");
        MethodVisitor constructorVisitor = cw.visitMethod(ACC_PUBLIC, "<init>", Type.getMethodDescriptor(Type.VOID_TYPE, new Type[]{}), null, null);
        constructorVisitor.visitCode();
        constructorVisitor.visitVarInsn(ALOAD, 0);
        constructorVisitor.visitMethodInsn(INVOKESPECIAL, superClassName, "<init>", Type.getMethodDescriptor(Type.VOID_TYPE, new Type[]{}));
        constructorVisitor.visitInsn(RETURN);
        constructorVisitor.visitMaxs(1, 1);
        constructorVisitor.visitEnd();

        //fields
        for (ClothoField field : fields) {
            log.trace("{}.visitField({},{},{},{},{})", cw, accessToOpcode(field.getAccess()), field.getName(), Type.getType(field.getType()).getDescriptor(), null, null);
            FieldVisitor fv = cw.visitField(accessToOpcode(field.getAccess()), field.getName(), Type.getType(field.getType()).getDescriptor(), null, null);
            //TODO: annotating embed vs reference
            if (field.isReference()) {
                log.trace("{}.visitAnnotation({},{}).visitEnd()", Type.getType(Reference.class).getInternalName(), true);
                fv.visitAnnotation(Type.getType(Reference.class).getInternalName(), true).visitEnd();
            }

            if (field.getConstraints() != null) {
                for (Constraint constraint : field.getConstraints()) {
                    String descriptor = Type.getDescriptor(constraint.constraintType);
                    Map<String, Object> values = constraint.values;
                    log.trace("{}.visitAnnotation({},{})", fv, descriptor, true);
                    AnnotationVisitor av = fv.visitAnnotation(descriptor, true);
                    //XXX: invalid annotation values cause hibernate validator to crash
                    for (String valueName : values.keySet()) {
                        handleAnnotationValue(av, valueName, values.get(valueName));
                    }
                    log.trace("{}.visitEnd", av);
                    av.visitEnd();
                }
            }

            log.trace("{}.visitEnd", fv);
            fv.visitEnd();
            //getters and setters
            if (field.getAccess() != Access.PRIVATE) {
                if (field.getAccess() != Access.READONLY) {
                    //setter
                    MethodVisitor mv = cw.visitMethod(ACC_PUBLIC, field.getSetterName(), Type.getMethodDescriptor(Type.VOID_TYPE, new Type[]{Type.getType(String.class)}),
                            null, null);
                    mv.visitCode();
                    mv.visitVarInsn(ALOAD, 0);
                    mv.visitVarInsn(ALOAD, 1);
                    mv.visitFieldInsn(PUTFIELD, this.getInternalName(), field.getName(), Type.getDescriptor(field.getType()));
                    mv.visitInsn(RETURN);
                    mv.visitMaxs(2, 2);
                    mv.visitEnd();
                }
                // getter
                MethodVisitor mv = cw.visitMethod(ACC_PUBLIC, field.getGetterName(), Type.getMethodDescriptor(Type.getType(field.getType()), new Type[]{}),
                        null, null);
                mv.visitCode();
                mv.visitVarInsn(ALOAD, 0);
                mv.visitFieldInsn(GETFIELD, this.getInternalName(), field.getName(), Type.getDescriptor(field.getType()));
                mv.visitInsn(ARETURN);
                mv.visitMaxs(1, 1);
                mv.visitEnd();
            }
        }

        //TODO: methods


        /*for (Function method : methods){
         mv = cw.visitMethod(method.getAccess(), method.getName(), Type.getType(method.getReturnType()).getInternalName(), method.getSignature(), null)
         //if method name is isValid, attach validator annotation
         //also append all field-level validators
         //bind handle to function call
         }*/

        cw.visitEnd();
        return cwriter.toByteArray();
    }

    public static int accessToOpcode(Access access) {
        if (access == Access.PUBLIC) {
            return ACC_PUBLIC;
        } else {
            return ACC_PROTECTED;
        }
    }

    private void handleAnnotationArray(AnnotationVisitor arrayVisitor, Object[] a) {
        for (Object value : a) {
            handleAnnotationValue(arrayVisitor, null, value);
        }
        log.trace("{}.visitEnd()", arrayVisitor);
        arrayVisitor.visitEnd();
    }

    private void handleAnnotationValue(AnnotationVisitor visitor, String name, Object value) {
        if (value.getClass().isArray()) {
            log.trace("{}.visitArray({})", visitor, name);
            handleAnnotationArray(visitor.visitArray(name), (Object[]) value);
        } else if (value instanceof Enum) {
            log.trace("{}.visitEnum({},{},{})", visitor, Type.getDescriptor(value.getClass()), value.toString());
            visitor.visitEnum(name, Type.getDescriptor(value.getClass()), value.toString());
        } //XXX: not handling the annotation case right now; no use case at the moment
        else {
            log.trace("{}.visit({},{})", visitor, name, value);
            visitor.visit(name, value);
        }
    }
}
