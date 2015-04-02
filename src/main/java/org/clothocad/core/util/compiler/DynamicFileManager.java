/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util.compiler;

import java.io.IOException;
import javax.tools.FileObject;
import javax.tools.ForwardingJavaFileManager;
import javax.tools.JavaFileObject;
import javax.tools.JavaFileObject.Kind;
import javax.tools.StandardJavaFileManager;

/**
 * From http://www.javablogging.com/dynamic-in-memory-compilation/
 */
public class DynamicFileManager extends
		ForwardingJavaFileManager {
    /**
    * Instance of JavaClassObject that will store the
    * compiled bytecode of our class
    */
    private ByteFileObject jclassObject;

    /**
    * Will initialize the manager with the specified
    * standard java file manager
    *
    * @param standardManger
    */
    public DynamicFileManager(StandardJavaFileManager
        standardManager) {
        super(standardManager);
    }

    /**
    * Will be used by us to get the class loader for our
    * compiled class. It creates an anonymous class
    * extending the SecureClassLoader which uses the
    * byte code created by the compiler and stored in
    * the JavaClassObject, and returns the Class for it
    */
    @Override
    public ClassLoader getClassLoader(Location location) {
        //TODO: hook into dynamic classloader
        return super.getClassLoader(location);
    }

    /**
    * Gives the compiler an instance of the JavaClassObject
    * so that the compiler can write the byte code into it.
    */
    @Override
    public JavaFileObject getJavaFileForOutput(Location location,
                                               String className,
                                               Kind kind,
                                               FileObject sibling)
    throws IOException{
        jclassObject = new ByteFileObject(className, kind);
        return jclassObject;
    }
    
    public ByteFileObject getOutputFile(){
        return jclassObject;
    }
}
