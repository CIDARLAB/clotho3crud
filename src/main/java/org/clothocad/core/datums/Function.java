package org.clothocad.core.datums;

import org.clothocad.core.datums.util.Language;
import org.clothocad.core.execution.Script;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.model.Person;

import lombok.Getter;

import java.util.List;
import java.util.Map;

//variable number of type parameterization?
//java object for first-class functions?
//declaration of requrirements?
public class Function extends Module {

    public Function(){};

    public Function(String name, Argument[] arguments, Class output, String source, Language language) {
        this.outputType = output;
        this.setName(name);
        this.args = arguments;
        this.language = language;
        if (language != null) {
            this.setCode(source, language);
        }
    }

    @Reference
    private Person author;

    @Getter 
    private Argument[] args;

    @Getter 
    private FunctionTest[] tests;

    //XXX: losing some duck-typing style flexibility
    private Class outputType;
    //XXX: if single return type, could make things more typesafe java-side
    //XXX: I don't even know what to do with this
    //XXX: all our target languages have single return value, so multiple return value is undefined
        //XXX: python has automatic tuple unpacking, is that what is intended?

    //XXX: @Replace(encoder="encodePrecondition", decoder="decodePrecondition")
    private Script precondition;

    //TODO: convert to dict-style
    @Override
    public String getCodeToLoad() {
        return code.encapsulateModule(code.modularizeFunction(code.getSource()), getSetup());
    }


    public static class FunctionTest {
        private  List<Object> args;
        private  Object value;

        public FunctionTest(List<Object> argValues, Object expectedResult) {
            this.args = argValues; 
            this.value = expectedResult;
        }

        public FunctionTest(){};
    }

    public String encodePrecondition() {
        if (precondition == null) return null;
        return precondition.toString();
    }

    public void decodePrecondition(Map obj) {
        if (!obj.containsKey("precondition")) return;
        String source = obj.get("precondition").toString();
        precondition = Module.generateScript(source, Language.valueOf(obj.get("language").toString()));
    }
}
