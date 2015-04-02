/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.schema.Converters;
import org.clothocad.core.schema.Schema;

/**
 *
 * @author prashantvaidyanathan
 */
public class ConverterSchema extends SharableObjBase 
{
    public Schema convertTo;
    public Schema convertFrom;
    
}
