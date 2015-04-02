package org.clothocad.core.aspects.Interpreter;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.clothocad.core.persistence.Persistor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


class Learner{
    final static Logger logger = LoggerFactory.getLogger(Learner.class);
    public static Persistor persistor;

    /**
     * Provided the list of n-gram features and the action that they map to,
     * persist the pattern
     * 
     * @param features
     * @param action 
     */
    public static void learn(String cmd, Map<String, Object> action) {
        String[] features = Handler.convertToFeatures(Utilities.tokenize(cmd));
        for(String feature : features) {
            if (feature != null) {
                HashMap<String, Integer> indexEntry = getIndexEntry(feature);
                    if (indexEntry == null) {
                        indexEntry = new HashMap<String, Integer>();
                    }
                    int currvalue = 0;
                    if(indexEntry.containsKey(action)) {
                        currvalue = indexEntry.get(action);
                    }
                    currvalue++;
                    indexEntry.put(action.toString(), currvalue);
                    if(currvalue > 100) {
                        autocollapse(indexEntry);
                    }
                    saveIndexEntry(indexEntry, feature);
            }
        }
    }
    

    /**
     * Provided the list of n-gram features and the action that they map to,
     * persist the pattern
     * 
     * @param features
     * @param action 
     */
    public static void forget(String cmd, String action) {
    String[]features = Handler.convertToFeatures(Utilities.tokenize(cmd));
        for(String feature : features) {
            HashMap<String, Integer> indexEntry = getIndexEntry(feature);
            if (indexEntry != null) {
                if(!indexEntry.containsKey(action)) {
                    continue;
                }
                int currvalue = indexEntry.get(action);
                if(currvalue<=1) {
                    indexEntry.remove(action);
                    saveIndexEntry(indexEntry, feature);
                    continue;
                }
                currvalue--;
                indexEntry.put(action, currvalue);
                saveIndexEntry(indexEntry, feature);
            }
        }
    }
    
    /* Instantiating Learner's own scanner. This is temporary */
    private static StdIn inputReader = new StdIn();

    /**
     * Learns the user's choice.
     */
    public static String Inquire(TreeMap<String, Double> ActionList) {
        giveOptions(ActionList);
        int choice = inputReader.readInt();
        String chosenCmd;
        try {
            chosenCmd = (String) ActionList.keySet().toArray()[choice - 1];
        } catch (IndexOutOfBoundsException ie) {
            logger.error("Your choice was invalid: {}", ie);
            return null;
        }
        logger.info("Your choice was {}", chosenCmd);
        return chosenCmd;
    }

    /**
     * Prints out the list of options referenced by cardinal numbers.
     */
    private static void giveOptions(TreeMap<String, Double> ActionList) {
        int cardNum = 1;
        logger.info("Choose your best option");
        for (String action : ActionList.keySet()) {
            logger.info("{}. {}", cardNum, action);
            cardNum ++;
        }
    }

    /**
     * Convenience method to fetch the HashMap associated with this entry
     * @param feature 
     */
    public static HashMap<String, Integer> getIndexEntry(String feature) {
        HashMap<String, Integer> out = persistor.loadFeature(feature);
        return out;
    }

    /**
     * Convenience method to persist an IndexEntry
     * @param feature
     * @param HashMap 
     */
    private static void saveIndexEntry(HashMap<String, Integer> StoreGrams, String feature) {
        persistor.persistFeature(StoreGrams, feature);
    }

    /* TODO replace autocollapse's job with Decider. If decide to keep then replace 
       the many for loops with Iterator removes                                    */
    /**
     * Once the count for an action in an indexEntry exceeds 100 counts, this gets 
     * called. It normalized all the counts in the entry, so it won't change the 
     * answer, but it keeps the numbers low so that new information has a chance 
     to take hold in there and
     * so that excessively rare actions get thrown out to avoid bloat.  This helps keep the
     * size of the data here minimal.
     * 
     * @param indexEntry 
     */
    private static void autocollapse(HashMap<String, Integer> indexEntry) {
        Set<String> keys = indexEntry.keySet();
        
        //Calculate the total counts in this index map
        int total = 0;
        for(String str : keys) {
            total += indexEntry.get(str);
        }
        
        //Normalize everything so that the total is 100 and remove rare things
        HashSet<String> removable = new HashSet<String>();
        for(String str : keys) {
            double oldvalue = (double) indexEntry.get(str);
            double newvalue = Math.floor(50.0*oldvalue / total);
            if(newvalue<1) {
                removable.add(str);
            } else {
                int ival = (int) newvalue;
                indexEntry.put(str, ival);
            }
        }
        
        //Remove the removables
        for(String str : removable) {
            indexEntry.remove(str);
        }
    }
}
