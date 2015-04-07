package org.clothocad.core.datums.util;

import java.util.Arrays;
import java.util.List;

public enum Language {
    JAVASCRIPT("JavaScript", Arrays.asList(new String[]{"js"})),
    PYTHON("Python", Arrays.asList(new String[]{"py"})),
    JAVA("Java", Arrays.asList(new String[]{"java"})),
    JSONSCHEMA("JSON Schema", Arrays.asList(new String[]{}));

    private Language(final String s, final List<String> strList){
        prettyName = s;
        extensions = strList;
    }

    public List<String> extensions() {
        return extensions;
    }

    public final String prettyName;
    private final List<String> extensions;

}
