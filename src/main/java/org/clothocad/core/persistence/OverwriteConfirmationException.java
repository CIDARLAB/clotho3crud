/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import java.util.Set;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
public class OverwriteConfirmationException extends RuntimeException {
    public final Set<ObjBase> modifiedObjects;

    public OverwriteConfirmationException(Set<ObjBase> modifiedObjects) {
        this.modifiedObjects = modifiedObjects;
    }
    
}
