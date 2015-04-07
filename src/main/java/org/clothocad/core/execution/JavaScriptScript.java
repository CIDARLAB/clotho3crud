package org.clothocad.core.execution;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.IdUtils;

import lombok.Getter;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class JavaScriptScript implements Script {
    //needs arg names, function name

    //how well does one engine for all functions scale?
    //do we need to manage the number of global functions? 

    //TODO: each script function needs to execute in its own scope
    //but
    //script functions should be cached somehow

    public JavaScriptScript() {}

    public JavaScriptScript(String source) {
        this.source = source;
    }

    @Getter
    private String source;

    @Override
    public Set<ObjectId> findDependencies() {
       String strippedSource =  source.replaceAll("\\s+","");
       Pattern pattern = Pattern.compile( "clotho\\.load\\(\"([0-9a-f]*)\"\\);");
       Matcher matcher = pattern.matcher(source);
       Set<ObjectId> output = new HashSet<>();
       while (matcher.find()) {
           output.add(new ObjectId(matcher.group(1)));
       }
       return output;
    }

    @Override
    public String generateImports(Collection<ObjectId> imports) {
        StringBuilder builder = new StringBuilder();
        String format = "var %s = clotho.load(\"%s\");\n";
        for (ObjectId id : imports) {
            ObjBase obj = IdUtils.get(id);
            String name = obj.getName();
            builder.append(String.format(format, name, id.toString()));
        }
        builder.append("\n");
        return builder.toString();
    }

    @Override 
    public String toString() {
        return source;
    }

    @Override
    public String modularizeFunction(String code) {
        //semicolon check
        code = code.trim();
        if (!code.endsWith(";")) {
            code = code + ";";
        }
        String format = "(function() {"
                + "var f = %s"
                + "return f;"
                + "}());";
        return String.format(format, code);
    }

    @Override
    public String encapsulateModule(String code, String setupcode) {
        String format = "(function() {"
                + "%s"
                + "return %s"
                + "}());";
        return String.format(format, setupcode, code);
    }
}
