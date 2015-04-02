/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.List;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
@NoArgsConstructor
public abstract class DataSource extends ObjBase {
    //XXX: revisit architecture later
    
    public abstract List<? extends ObjBase> importData(String info);
    public abstract void exportData(List<ObjBase> exports);
    
    @Getter
    @Setter
    String url;
    
    protected boolean canImport = false;
    protected boolean canExport = false;
    
    public DataSource(String name, String url){
        super(name);
        this.url = url;
    }
            
}
