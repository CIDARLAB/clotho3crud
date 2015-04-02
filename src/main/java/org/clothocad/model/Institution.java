/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import javax.validation.constraints.NotNull;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.SharableObjBase;

/**
 *org.clothocad.model.Institution
 * @author jcanderson
 */
@NoArgsConstructor
public class Institution extends SharableObjBase {

    /**
     * Constructor from raw data
     * @param name
     * @param city
     * @param state
     * @param country
     */
    
    @Setter
    @Getter
    @NotNull
    private String city, state, country;
    
    //TODO:unique name constraint
    public Institution( String name, String city, String state, String country ) {
        super(name, null);
        this.city = city;
        this.state = state;
        this.country = country;
    }

}
