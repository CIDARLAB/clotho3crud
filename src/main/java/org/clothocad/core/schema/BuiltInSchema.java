/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.fasterxml.jackson.annotation.JsonView;
import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.util.HashSet;
import java.util.Set;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.jackson.JSONViews;

/**
 *
 * @author spaige
 */
@NoArgsConstructor
@Slf4j
@EqualsAndHashCode(callSuper = false, of = {"binaryName", "c"})
public class BuiltInSchema extends JavaSchema {

    @Reference
    protected BuiltInSchema superClass;
    public BuiltInSchema(Class<? extends ObjBase> c, BuiltInSchema superClass) {
        this(c);
        this.superClass = superClass;
    }

    public BuiltInSchema(Class<? extends ObjBase> c) {
        super();
        fields = extractFields(c);
        //TODO: functions

        binaryName = c.getCanonicalName();
        setName(c.getSimpleName());
        this.c = c;
        setId(new ObjectId(c.getCanonicalName()));
    }
    Class<? extends ObjBase> c;

    @Override
    public void setSource(String source) {
        throw new UnsupportedOperationException("You cannot set the source of a built-in schema.");
    }
    private String binaryName;

    @Override
    public String getBinaryName() {
        return binaryName;
    }

    @Override
    public Class<? extends ObjBase> getEnclosedClass(ClassLoader cl) {
        return c;
    }

    private static boolean fieldHasAnnotation(Field field, Class<? extends Annotation> annotation) {
        Annotation[] annotations = field.getAnnotations();
        for (int i = 0; i < annotations.length; i++) {
            if (annotations[i].annotationType().equals(annotation)) {
                return true;
            }
        }
        return false;
    }

    private static boolean fieldIsDBOnly(Field field) {
        Annotation[] annotations = field.getAnnotations();
        for (int i = 0; i < annotations.length; i++) {
            if (annotations[i].annotationType().equals(JsonView.class)) {
                Class[] views = ((JsonView) annotations[i]).value();
                return views.length == 1 && views[i].equals(JSONViews.Internal.class);
            }
        }
        return false;
    }

    private Set<ClothoField> extractFields(Class c) {
        Set<ClothoField> fields = new HashSet<>();
        for (Field field : c.getDeclaredFields()) {
            if (fieldIsDBOnly(field)) {
                continue;
            }
            //TODO: do not store transient fields
            fields.add(new ClothoField(field));
        }
        return fields;
    }
}
