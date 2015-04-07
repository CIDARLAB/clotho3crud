package org.clothocad.core;

import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.nosecurity.NoSecurityModule;
import org.clothocad.model.FreeForm;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;
import org.clothocad.model.SimpleSequence;

public class PersistorModularityDemo {
    public static void main(String[] args) {
        //Setup: make a builder configured with the modules you want, then get 
        //the classes you want to use
        ClothoBuilder builder = new ClothoBuilder(
                //We have to have a security manager of some kind - this module 
                //sets up one that just does 'pass-through' checks on permissions
                new NoSecurityModule(),
                //Sets up the persistor to use Jongo/MongoDB as the database
                new JongoModule());

        //get a persistor
        Persistor p = builder.get(Persistor.class);
        
        //use the persistor to do stuff
        Person demoPerson = new Person("Demo Person");
        Sequence demoSeq = new SimpleSequence("ATCG", demoPerson);
        Part part = new Part("Demo Part", "This is a demo part.", demoSeq, demoPerson);
        part.setFormat(new FreeForm());

        p.save(part);

        Part retrievedPart = p.get(Part.class, part.getId());
    }
}
