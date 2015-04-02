/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.execution;

/**
 *
 * @author prashantvaidyanathan
 */
public class SequenceConverters {
    
    public static String getSequence(String sequence)
    {
        return sequence;
    }
    public static String revString(String sequence)
    {
        String res = new StringBuilder(sequence).reverse().toString();
        return res;
    }
    
}
