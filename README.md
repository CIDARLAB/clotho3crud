Clotho 
===========================================
CRUD (Create, Read, Update, Delete)
Clotho CRUD is a Java library that helps manage the data used for building synthetic biology software.

Data in Clotho must adhere to specified Schemas, which may represent models from Biological Parts to entire Experimental Designs. This standardized method of storing and handling data allows for easy usage in your programs, letting you focus on everything else.

Installation
------------
Simply download the latest jar file in /distributions and add it to your project. The jar contains all the necessary dependencies for Clotho CRUD.

Requirements
------------
[mongoDB](https://www.mongodb.org/downloads) is required. Currently versions up to 2.4 are supported.

Quick Start
-----------
Step 1: Import the relevant files into your Java class
```
// Required in order to build Clotho
import org.clothocad.core.ClothoBuilder;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.nosecurity.NoSecurityModule;

// Schemas needed for our demo example
import org.clothocad.model.FreeForm;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;
import org.clothocad.model.SimpleSequence;
```
Step 2: Build Clotho with the modules you want and get a Persistor
```
ClothoBuilder builder = new ClothoBuilder(
        // We have to have a security manager of some kind - this module 
        // sets up one that just does 'pass-through' checks on permissions
        new NoSecurityModule(),
        // Sets up the persistor to use Jongo/MongoDB as the database
        new JongoModule());
// Persistor will be used for CRUD operations
Persistor p = builder.get(Persistor.class);
```
Step 3: Create, Read, Update, and Delete!
```
Person demoPerson = new Person("Demo Person");
Sequence demoSeq = new SimpleSequence("ATCG", demoPerson);
Part part = new Part("Demo Part", "This is a demo part.", demoSeq, demoPerson);
part.setFormat(new FreeForm());

// Save the created Schema to the database
p.save(part);

// Get Schemas from the database
Part retrievedPartById = (Part) p.get(part.getId());

// Delete Schema from the database
p.delete(part.getId());
```

Documentation
------------
In progress (see Wiki)

Contact
-------
If you want to be informed about new releases, bug fixes, and general
updates about Clotho, subscribe to the [Clotho Users Group]( 
https://groups.google.com/group/clotho-users)

Contributing to this project
----------------------------
Anyone and everyone is welcome to contribute. Please take a moment to
review the [guidelines for contributing](CONTRIBUTING.md).

* [Bug reports](CONTRIBUTING.md#bugs)
* [Feature requests](CONTRIBUTING.md#features)
* [Pull requests](CONTRIBUTING.md#pull-requests)

License
-------
Please see [license.txt](/license.txt)
