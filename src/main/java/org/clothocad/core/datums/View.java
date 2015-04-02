package org.clothocad.core.datums;

import java.util.List;
import java.util.Map;
import javax.validation.constraints.AssertTrue;

import org.clothocad.core.datums.util.ClothoField;

import org.clothocad.model.Person;


/**
 * A View is the factory for HTML widgets
 * 
 * It is one of the 4 primary Sharables
 * 
 * @author John Christopher Anderson
 */
public class View extends SharableObjBase {

    private View(Person author, 
            String name, 
            String description,
            List<ClothoField> inputArguments,
            Function canUpdate, 
            String graphicsScript,
            String onShowScript,
            String onUpdateScript) {
    	
    	super(name, author);
        
        this.inputArguments = inputArguments;
        this.canUpdate = canUpdate;
        this.graphicsScript = graphicsScript;
        this.onShowScript = onShowScript;
        this.onUpdateScript = onUpdateScript;
        this.name = name;
        this.description = description;
        
    }

    public static View create(
            Person author, 
            String name, 
            String description,
            List<ClothoField> inputArgs,
            Function canUpdate, 
            String graphicsScript,
            String onShowScript,
            String onUpdateScript) {
        
    	View view = new View(author, name, description, 
    			inputArgs, canUpdate, graphicsScript, 
    			onShowScript, onUpdateScript);
    	
        //Collector.get().add(view);
        return view;
    }
    
    /***
    @Override
    public Map<String, Object> toJSON() {
        try {
            JSONSerializer serializer = new JSONSerializer().exclude("*.class");
            serializer.prettyPrint(true);
            String serial = serializer.deepSerialize( this );
            return new Map<String, Object>(serial);
        } catch (JSONException ex) {
            return null;
        }
    }
	***/
    


    public int getInstanceCount() {
        return instanceCount;
    }

    public String getLargeIconURL() {
        return largeIconURL;
    }

    public String getSmallIconURL() {
        return smallIconURL;
    }

    public List<ClothoField> getInputArguments() {
        return inputArguments;
    }

    public Function getCanUpdate() {
        return canUpdate;
    }

    public String getOnUpdateScript() {
        return onUpdateScript;
    }
    
    public String getGraphicsScript() {
        return this.graphicsScript;
    }
    
    public String getOnShowScript() {
        return this.onShowScript;
    }

    //Still need to implement this:
    private List<ClothoField> inputArguments;
    private Function canUpdate;
    private String onUpdateScript;
    private String graphicsScript;
    private String onShowScript;
    

    //Metadata
    private String id;
    private String name;
    private String description;
    
    private String authorId;
    private String smallIconURL;
    private String largeIconURL;
    
    private int instanceCount = 0;

    @AssertTrue
    public boolean validate(Map<String, Object> obj) {
        return true;
    }
}
