package org.clothocad.model;

import org.clothocad.core.datums.ObjBase;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.List;

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
