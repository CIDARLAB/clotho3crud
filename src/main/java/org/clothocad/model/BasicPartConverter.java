package org.clothocad.model;

import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.InferredSchema;
import org.clothocad.core.schema.Schema;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author spaige
 */
public class BasicPartConverter extends Converter<Part> {
    static Set<String> names = new HashSet<String>();
    static {
        names.add("eugene.dom.components.Part");
    }

    public BasicPartConverter(Persistor p) {
        super(p.get(Schema.class, new ObjectId("org.clothocad.model.Part")), new HashSet<Schema>(), names);
    }

    @Override
    protected Part guardedConvert(Map data, String schemaName) {
        switch (schemaName){
            case "eugene.dom.components.Part":
                return convertEugenePartToBasicPart(data);
            default:
                return null;
        }
        
    }

    public static Part convertEugenePartToBasicPart(Map<String, Object> eugenePart) {
        Person author = new Person("Anonymous");
        Sequence partSeq = new SimpleSequence(eugenePart.get("Sequence").toString(), author);
        Part part = new Part(eugenePart.get("Name").toString(), partSeq, author);
        if (eugenePart.containsKey("_id")) 
            part.setId(new ObjectId(eugenePart.get("_id").toString()));
        try {
            Feature feature = new Feature(eugenePart.get("Name").toString(), 
                    Feature.FeatureRole.valueOf(eugenePart.get("PartType").toString().toUpperCase()),
                    author);
            feature.setSequence(partSeq);
            Annotation seqAnnotation = partSeq.createAnnotation(eugenePart.get("Name").toString(),
                    0, partSeq.getSequence().length() - 1, true, author);
            seqAnnotation.setFeature(feature);
        } catch (IllegalArgumentException e) {
        }
        return part;
    }

    @Override
    protected Part guardedConvert(Map data, Schema type) {
        if (type instanceof InferredSchema ){
            return guardedConvert(data, type.getName());
        }
        
        throw new UnsupportedOperationException("No class schemas supported"); //To change body of generated methods, choose Tools | Templates.
    }
}
