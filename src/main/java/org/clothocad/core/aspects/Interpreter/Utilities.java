package org.clothocad.core.aspects.Interpreter;

import java.util.Comparator;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

public class Utilities {
    /**
     * Recursive factorial function.
     * @param num to calculate it's factorial
     * @return num! integer factorial of the input
     */
    public static int fact(int num) {
        if (num <= 0) {
            return 1;
        } else {
            return num * fact(num - 1);
        }
    }

    /**
     * Combinations
     * @param n, r
     * @return num of combinations
     */
    public static int comb(int num, int r) {
        return fact(num)/(fact(r) * fact(num-r));
    }

    /**
     * Tokenize the raw cmd into a String array of ordered words, then
     * run any substitution of proper names
     * @param cmd
     * @return 
     */
    public static String[] tokenize(String cmd) {
        //Break the cmd into words
        String[] out = cmd.toLowerCase().split("[\\s,.;]+"); 
        
        //Exchange out proper names for Schema uuids via Collator
        for(int i = 0; i < out.length; i++) {
            String newword = out[i];
            
            /* If the word was a proper name, replace it with 
             * the Schema or "assistant" reference
             */
            if(!out[i].equals(newword)) {
                out[i] = newword;
            }   
        }
        return out;
    }
}
