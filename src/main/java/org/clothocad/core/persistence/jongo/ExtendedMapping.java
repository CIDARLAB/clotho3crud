/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectReader;
import org.jongo.marshall.jackson.configuration.Mapping;
import org.jongo.marshall.jackson.configuration.ReaderCallback;
import org.jongo.marshall.jackson.configuration.WriterCallback;

/**
 *
 * @author spaige
 */
class ExtendedMapping extends Mapping{

    public ExtendedMapping(ObjectMapper mapper, ReaderCallback readerCallback, WriterCallback writerCallback) {
        super(mapper, readerCallback, writerCallback);
    }

    ObjectReader getUpdatingReader(Class clazz) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
