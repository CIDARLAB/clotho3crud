/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.communication.imports.IceImporter;

/**
 *
 * @author spaige
 */
@NoArgsConstructor
public class ICEDataSource extends DataSource{
    
    
    public ICEDataSource(String name, String url){
        super(name, url);
    }

    public boolean canImport = true;
    
    @Override
    public List<? extends ObjBase> importData(String info) {
        String[] idStrings = info.split(",");
        List<Integer> ids = new ArrayList<>();
        for (String idString : Arrays.asList(idStrings) ){
            ids.add(Integer.parseInt(idString.trim()));
        }
        return IceImporter.importData(ids);
    }

    @Override
    public void exportData(List<ObjBase> exports) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
