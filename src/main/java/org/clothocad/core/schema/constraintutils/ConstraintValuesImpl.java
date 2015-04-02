/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema.constraintutils;

import java.util.HashSet;
import java.util.Set;
import javax.validation.Payload;

/**
 *
 * @author spaige
 */
public class ConstraintValuesImpl implements ConstraintValues {
    Class<?>[] groups;
    String message;
    Class<? extends Payload> payload;

    @Override
    public Object get(String key) {
        switch(key){
            case("groups"):
                return groups;
            case("message"):
                return message;
            case("payload"):
                return payload;
            default:
                throw new IllegalArgumentException();
        }
    }

    @Override
    public void put(String key, Object value) {
        switch(key){
            case("groups"):
                groups = (Class<?>[]) value;
                break;
            case("message"):
                message = (String) value;
                break;
            case("payload"):
                payload = (Class<? extends Payload>) payload;
                break;
            default:
                throw new IllegalArgumentException();
        }
    }

    @Override
    public Set<String> keySet() {
        Set<String> out = new HashSet<>();
        if (message != null) out.add("message");
        if (payload != null) out.add("payload");
        if (groups != null) out.add("groups");
        return out;
    }
}
