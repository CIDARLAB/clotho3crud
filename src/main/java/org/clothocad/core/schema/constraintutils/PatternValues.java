/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema.constraintutils;

import java.util.List;
import java.util.Set;
import javax.validation.constraints.Pattern;

/**
 *
 * @author spaige
 */
public class PatternValues extends ConstraintValuesImpl{
    private Pattern.Flag[] flags;
    private String regexp;

    @Override
    public Object get(String key) {
        switch(key){
            case("flags"):
                return flags;
            case("regexp"):
                return regexp;
        }
        return super.get(key);
    }

    @Override
    public void put(String key, Object value) {
        switch(key){
            case("flags"):
                if (value instanceof List) {
                    flags = new Pattern.Flag[((List) value).size()];
                    for (int i = 0; i < flags.length; i ++){
                        flags[i] = Pattern.Flag.valueOf(((List) value).get(i).toString());
                    }
                }
                else flags = (Pattern.Flag[]) value;
                break;
            case("regexp"):
                regexp = (String) value;
                break;
            default:
                super.put(key, value);
        }

    }

    @Override
    public Set<String> keySet() {
        Set<String> out = super.keySet();
        if (flags != null) out.add("flags");
        if (regexp != null) out.add("regexp");
        return out;
    }
}
