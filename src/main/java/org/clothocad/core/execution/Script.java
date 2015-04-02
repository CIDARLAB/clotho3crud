/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import java.util.Collection;
import java.util.Set;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public interface Script{
    public Set<ObjectId> findDependencies();
    
    public String getSource();

    public String generateImports(Collection<ObjectId> listedButNotDeclared);
    
    public String modularizeFunction(String code);
    
    public String encapsulateModule(String code, String setupCode);
}
