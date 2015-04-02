/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */

package org.clothocad.core.schema;

import static java.util.Arrays.asList;
import javax.tools.JavaCompiler;
import javax.tools.JavaCompiler.CompilationTask;
import javax.tools.ToolProvider;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.util.compiler.DynamicFileManager;
import org.clothocad.core.util.compiler.JavaSourceFromString;
import org.clothocad.model.Person;
import org.objectweb.asm.ClassReader;


public class JavaSchema 
		extends Schema {

    public JavaSchema() {super();}
	
    public JavaSchema(Person author,
            String name, 
            String description,
            String source) {
    	
    	super(name, description, author);
        
        setSource(source);
    } 
    
    
    protected byte[][] innerClasses;
    
    //TODO:
    public JavaSchema(String source){
        this.setSource(source);
        this.extractMetadata();
    }
    
    //TODO: handle files that result in multiple source files;
    //TODO: supply uuid's of referenced classes to resolve name conflicts
    @Override
    public void setSource(String source){
        this.source = source;
        
        JavaSourceFromString sourceFile = new JavaSourceFromString(getName(), source);
        JavaCompiler compiler = ToolProvider.getSystemJavaCompiler();
        DynamicFileManager fileManager = new DynamicFileManager(compiler.getStandardFileManager(null, null, null));
        CompilationTask task = compiler.getTask(null, fileManager, null, null, null, asList(sourceFile));

        boolean status = task.call();

        if (!status){
            //error handling
            return;
        } 
        
        this.classData = fileManager.getOutputFile().getBytes();
        //patch references to other db-residing classes to uuid instead of 'pretty' id
        //TODO: how to manage external dependencies?
        //complain if non objbase class

    }
    
    private void extractMetadata(){
        ClassReader reader = new ClassReader(classData);
        reader.accept(new ClassParser(this), 0);   
        //fields
        //methods
        //author/name/description 
    }

    @Override
    public Language getLanguage() {
        return Language.JAVA;
    }
    
}

