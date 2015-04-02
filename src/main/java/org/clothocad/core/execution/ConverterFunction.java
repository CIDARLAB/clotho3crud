/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.execution;

import org.clothocad.core.datums.Argument;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.schema.Schema;
import org.clothocad.model.ConverterSchema;

/**
 *
 * @author prashantvaidyanathan
 */
public class ConverterFunction extends ConverterSchema
{
    Function convFunction;
    ConverterFunction()
    {
        this.convFunction = new Function();
        
    }
    ConverterFunction(Schema convTo, Schema convFrom, Function func)
    {
        this.convertFrom = convFrom;
        this.convertTo = convTo;
        this.convFunction = func;
    }
    public Schema getConvertToSchema()
    {
            return convertTo;
    }
    public Schema getConvertFromSchema()
    {
            return convertFrom;
    }
    public Function getFunction()
    {
        return convFunction;
    }
    
    
}
