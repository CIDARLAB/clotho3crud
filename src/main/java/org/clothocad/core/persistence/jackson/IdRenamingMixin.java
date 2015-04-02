/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.annotation.JsonProperty;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public abstract class IdRenamingMixin {
    @JsonProperty("_id") private ObjectId id;
}
