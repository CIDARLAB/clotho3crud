package org.clothocad.core.aspects.Interpreter;

import java.util.Comparator;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Iterator;

import static org.clothocad.core.aspects.Interpreter.Utilities.comb;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

class Handler {
    public static Logger logger = LoggerFactory.getLogger(Handler.class);
     /**
     * Fetch the list of actions with unsorted score and return them in a
     * TreeMap sorted by Value. All of the values must sum up to 1.
     * @param features
     * @return 
     */
    public static TreeMap<String, Double> query(String[] features) {
        TreeMap<String, Double> populator = new TreeMap<String, Double>();
        /* For each feature of this command, normalize and gather up the scores. */
        for(String feature : features) {
            HashMap<String, Integer> StoredGrams = Learner.getIndexEntry(feature);
            if (StoredGrams != null) { 
                populatingActionsList(StoredGrams, populator);
            }
        }
        /* Normalize and put into the value-sorted TreeMap */
        ValueComparator comparator = new ValueComparator(populator);
        TreeMap<String, Double> sortedOut = new TreeMap<String, Double>(comparator);
        double totalScore = 0;
        for (double score : populator.values()) {
            totalScore += score;
        }
        for (String action : populator.keySet()) {
            double normScore = populator.get(action) / totalScore;
            sortedOut.put(action, normScore);
        }
        checkProb(sortedOut);
        return Threshold(sortedOut);
    }

    /**
     * Populates out with StoredGram scores normalizing them previous to
     * putting it in.
     * @params HashMap StoredGrams, TreeMap out
     */
     private static void populatingActionsList(HashMap<String, Integer> StoredGrams,
                                               TreeMap<String, Double> out) {
        double totalScore = 0;
        for(String action : StoredGrams.keySet()) {
            totalScore += StoredGrams.get(action);
            if (out.containsKey(action)) {
                totalScore += out.get(action);
            }
        }
        for(String action : StoredGrams.keySet()) {
            double newScore = StoredGrams.get(action);
            if (out.containsKey(action)) {
                newScore += out.get(action);
            }
            out.put(action, newScore/totalScore);
        }
     }

    /**
     * Provided a list of words, it computes the combinations of 
     * pair and triple of words that will be the features upon 
     * which the command gets indexed.
     * @param words
     * @return 
     */
    public static String[] convertToFeatures(String[] sentence) {
        /* Generating array for features */
        int n = sentence.length;
        int pairCombs = comb(n, 2);
        int tripleCombs = comb(n, 3);
        int totalCombs = pairCombs + tripleCombs;

        /* +1 is for the wholeword */
        String[] out = new String [totalCombs + 1];
        
        int index = 0;

        /* Placing the pairs */
        if (n > 2) {
            for(int i = 0; i < (sentence.length - 1); i++) {
                for(int j = i+1; j < sentence.length; j++) {
                    String joinWord = sentence[i];
                    /* TODO testing this taken out
                    if(j-i == 1) {
                        joinWord += ".";
                    } else {
                        joinWord += "_";
                    }*/
                    joinWord += ".";
                    joinWord += sentence[j];
                    out[index] = joinWord;
                    index += 1;
                }
            }
        }

        if (n > 3) {
            /* Placing the triples */
            for(int i = 0; i < (sentence.length - 2); i++) {
                for(int j = i+1; j < (sentence.length - 1); j++) {
                    for(int k = j+1; k < sentence.length; k++) {
                        String joinWord = sentence[i];
                        if(j-i == 1) {
                            joinWord += ".";
                        } else {
                            joinWord += "_";
                        }
                        joinWord += sentence[j];
                        
                        if(k-j == 1) {
                            joinWord += ".";
                        } else {
                            joinWord += "_";
                        }
                        joinWord += sentence[k];
                        out[index] = joinWord;
                        index += 1;
                    }
                }
            }
        }
        /* Placing the whole sentence */
        String wholeword = sentence[0];
        for(int i=1; i<sentence.length; i++) {
            wholeword+=".";
            wholeword+=sentence[i];
        }
        out[index] = wholeword;
        return out;    
    }

    /**
     * Invariant to check that probs are calculated correctly
     */
    private static void checkProb(TreeMap<String, Double> ActionList) {
        double total = 0;
        for (Map.Entry<String, Double> entry : ActionList.entrySet()) {
            total += entry.getValue();
        }
        try {
            if (total != 0 && total != 1) {
            throw new Exception();
            }
        } catch (Exception e) {
            logger.error(
                       "Probability badly calculated, values summed up to" + total,
                       e);
        }
    }

    /**
     * Keeps the more possible options and trims the rare ones. For now
     * the threshold is set to everything within 0.4.
     */
    private static TreeMap<String,Double> Threshold(TreeMap<String, Double> ActionList) {
        if (ActionList.firstEntry() != null) {
            double current = ActionList.firstEntry().getValue();
            double previous = current;
            double threshold = 0.4;
            for (Iterator<Map.Entry<String, Double>> itr = ActionList.entrySet().iterator();
                    itr.hasNext();) {
                Map.Entry<String, Double> entry = itr.next();
                current = entry.getValue();
                if ((previous - current) > threshold) {
                    while (itr.hasNext()) {
                        itr.remove();
                        itr.next();
                    }
                    itr.remove();
                }
                previous = current;
            }
        }
        return ActionList;
    }
}

/**
 * Comparator used to sort a Map by its values
 */
class ValueComparator implements Comparator<String> {
    Map <String, Double> base;
    public ValueComparator(Map <String, Double> base){
        this.base = base;
    }
    
    public int compare(String a, String b) {
        if (base.get(a) < base.get(b)) {
            return 1;
        } else if (base.get(a) == base.get(b)) {
            return 0;
        } else {
            return -1;
        }
    }
}
