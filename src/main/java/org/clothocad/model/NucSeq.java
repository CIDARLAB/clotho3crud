/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import java.io.Serializable;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.validation.constraints.NotNull;
import lombok.Getter;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;

@NoArgsConstructor
public class NucSeq 
		extends ObjBase 
		implements Serializable {

    /**
     * Constructor for a NucSeq (dna sequence) object.  The sequence ends up in the
     * private variable 'theSequence'.  In constructing it, it checks whether it is
     * a legitimate RNA or DNA molecule and whether it contains degenerate codons,
     * which is stored as a boolean 'isDegenerate'.
     *
     * @param inputSeq
     * @param strandedness
     * @param circularity
     */
    
    private boolean isSingleStranded, isCircular;
    //are linear and circular mutually exclusive?
    private boolean isDegenerate, isLinear, isRNA;
    private boolean isLocked;
    
    @Getter
    private Set<Annotation> annotations;
    
    @NotNull
    private String sequence;
    
    public NucSeq( String inputSeq, boolean strandedness, boolean circularity ) {
        super("nucseq");

        isSingleStranded = strandedness;
        isCircular = circularity;
        lowerArray = new boolean[ inputSeq.length() ];
        if ( !initiateNucSeq( inputSeq ) ) {
            return;
        }


        //NEED TO USE A DIFFERENT OBJBASE CONSTRUCTOR TO SET HASH AS UUID
        //_myUUID = generateUUIDAsHash(getSeq());
    }

    //alternate constuctor if circularity and strandedness isn't specified
    public NucSeq( String inputSeq ) {
        this( inputSeq, false, false );
    }
    
    private static final ImmutableList<String> START_CODONS = ImmutableList.of(
            "ATG", 
            "GTG",
            "TGG",
            "RTG");
    
    private static final ImmutableList<String> STOP_CODONS = ImmutableList.of(
            "TAA",
            "TAG",
            "TGA",
            "TRA",
            "TAR");
        
   

    public boolean initiateNucSeq( String inputSeq ) {
        
        if(inputSeq == null) {
            sequence = null;
            return false;
        }

        char currentchar;
        StringBuffer seq = new StringBuffer();
        lowerArray = new boolean[ inputSeq.length() ];

        //Check whether this is an RNA, DNA, or a bad seq, and if is degenerate:
        loopy:
        for ( int i = 0; i < inputSeq.length(); i++ ) {
            currentchar = inputSeq.charAt( i );
            char upperChar = Character.toUpperCase( currentchar );

            //Put the case in a format array
            if ( currentchar == upperChar ) {
                lowerArray[i] = false;
            } else {
                lowerArray[i] = true;
            }

            //Build up the new format-free sequence
            seq.append( upperChar );

            switcheroo:
            switch ( upperChar ) {
                case '.':
                    isLinear = true;
                    break;
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                    break;
                case 'B':
                case 'D':
                case 'H':
                case 'K':
                case 'M':
                case 'N':
                case 'R':
                case 'S':
                case 'V':
                case 'W':
                case 'Y':
                    isDegenerate = true;
                    break;
                case 'U':
                    isRNA = true;
                    break;
                default:
                    System.out.println( "Nucseq had to break on " + upperChar );
                    sequence = null;
                    return false;
            }
        }
        sequence = seq.toString();
        return true;
    }
    
    
    public Map<String,Object> getNucSeqMap() throws IllegalArgumentException, IllegalAccessException
    {
        Map<String,Object> nucsmap = new HashMap<String,Object>();
        
        Field nseqFields[] = NucSeq.class.getDeclaredFields();
        for(Field nsfield:nseqFields)
        {
            nsfield.setAccessible(true);
            nucsmap.put(nsfield.getName(), nsfield.get(this));
        }
        
        return nucsmap;
    }

     /**
     * Finds indices of Open Reading Frames in a given nucleotide sequence and
     * returns them as HashMap with start indices for keys and end indices for
     * values.
     *
     * @param s String to check for ORFs
     * @param forward Boolean, set to 'true' for finding forward reading frames
     *          or to 'false' to find ORFs in the reverse complement
     */
    @SuppressWarnings (value="unchecked")
    public HashMap<Integer, Integer> findORFs(boolean forward, boolean multipleStartCodons) {
        String seq = sequence;
        int len = seqLength();
        HashMap orfs = new HashMap();
        if (isCircular()) {
            seq = seq.concat(seq);
        }

        Pattern pattern = Pattern.compile(makeORFRegExp(multipleStartCodons, isDegenerate()), Pattern.CASE_INSENSITIVE);
        Matcher matcher;
        if (forward) {
            matcher = pattern.matcher(seq);
        }
        else {
            matcher = pattern.matcher(revComp());
        }

        int end;
        int start;
        int pos = 0;
        while (matcher.find(pos)) {
            start = matcher.start();
            end = matcher.end();

            if (end > len) {
                end = end - len;
            }

            if (!(start >= len || matcher.group().length() > len)) {
                if (forward) {
                    orfs.put(start, end);
                }
                else {
                    orfs.put(len - start, len - end);
                }
            }

            pos = matcher.start() + 3;
        }

        return orfs;
    }

        /**
     * Returns a regular expression for finding Open Reading Frames using
     * data from the codon table
     */
    private String makeORFRegExp(boolean multipleStartCodons, boolean allowDegen) {
        String regexp = "(";
        if (multipleStartCodons) {
            for (int i = 0; i < START_CODONS.size(); i++) {
                if (i + 1 < START_CODONS.size()) {
                    regexp = regexp + seqToRegExp(START_CODONS.get(i), allowDegen) + "|";
                }
                else {
                    regexp = regexp + seqToRegExp(START_CODONS.get(i), allowDegen) + ")";
                }
            }
        }
        else {
            regexp = regexp + seqToRegExp(START_CODONS.get(0), allowDegen) + ")";
        }
        regexp = regexp + "(...)*?(";
        for (int i = 0; i < STOP_CODONS.size(); i++) {
            if (i + 1 < STOP_CODONS.size()) {
                regexp = regexp + seqToRegExp(STOP_CODONS.get(i), allowDegen) + "|";
            }
            else {
                regexp = regexp + seqToRegExp(STOP_CODONS.get(i), allowDegen) + ")";
            }
        }
        return regexp;
    }

    
    //TODO: cache regexified versions of start and stop codons
    /**
     * Takes a sequence and transforms it into a regular expression for
     * searches.
     *
     * @param s String containing the sequence
     */
    private String seqToRegExp(String seq, boolean degen) {
        String regexp = "";
        String c;
        String prefix = "";
        String suffix = "";
        int rootstart = 0;
        int rootend = seq.length();
        if (seq.indexOf("<") != -1   ||   seq.indexOf(">") != -1) {
            for (int i = 0; i < seq.length(); i++) {
                c = seq.substring(i,i+1);
                if (c.equalsIgnoreCase("<")) {
                    Pattern pattern = Pattern.compile(".", Pattern.CASE_INSENSITIVE);
                    Matcher matcher = pattern.matcher(regexp);
                    while (matcher.find()) {
                        prefix = "(" + prefix + matcher.group() + ")?";
                    }
                    rootstart = i+1;
                }
                else if (c.equalsIgnoreCase(">")) {
                    Pattern pattern = Pattern.compile(".", Pattern.CASE_INSENSITIVE);
                    Matcher matcher = pattern.matcher(seq.replaceFirst(regexp,""));
                    int counter = 0;
                    matcher.find();
                    while (matcher.find()) {
                        suffix = suffix + "(" + matcher.group();
                        counter++;
                    }
                    for (int j=0; j<counter; j++) {
                        suffix = suffix + ")?";
                    }
                    rootend = i;
                }
                else {
                    regexp = regexp + c;
                }
            }
        }

        //System.out.println("Is");

        regexp = prefix + seq.substring(rootstart, rootend) + suffix;

        if (degen) {
            //The following block strings the regex and protects with @ symbol
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[aA][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[aA]", "[A@D@H@M@N@R@V@W]");
            //System.out.println("here?C " + regexp.length());
            // (regexp.matches("[a-z[A-Z]\\[\\]]*\\]?(?<!@)[cC]\\[?[a-z[A-Z]\\[\\]]*"))
                regexp = regexp.replaceAll("(?<!@)[cC]", "[@BC@H@M@N@S@V@Y]");
            //System.out.println("here?G " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]\\]]*(?<!@)[gG][a-z[A-Z]\\[]*"))
                regexp = regexp.replaceAll("(?<!@)[gG]", "[@B@DG@K@N@R@S@V]");
            //System.out.println("here?T " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]\\]]*(?<!@)[tT][a-z[A-Z]\\[]*"))
                regexp = regexp.replaceAll("(?<!@)[tT]", "[@B@D@H@K@NT@W@Y]");
            //System.out.println("here?U " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[uU][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[uU]", "[@B@D@H@K@NU@W@Y]");
            //System.out.println("here?B " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[bB][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[bB]", "[BC@DG@H@K@M@N@R@STU@V@W@Y]");
            //System.out.println("here?D " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[dD][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[dD]", "[ABDG@H@K@M@N@R@STU@V@W@Y]");
            //System.out.println("here?H " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[hH][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[hH]", "[ABCDH@K@M@N@R@STU@V@W@Y]");
            //System.out.println("here?K " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[kK][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[kK]", "[BDGHK@N@R@STU@V@W@Y]");
            //System.out.println("here?M " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[mM][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[mM]", "[ABCDHM@N@R@S@V@W@Y]");
            //System.out.println("here?N " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[nN][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[nN]", "[ABCDGHKMN@R@STU@V@W@Y]");
            //System.out.println("here?R " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[rR][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[rR]", "[ABDGHKMNR@S@V@W]");
            //System.out.println("here?S " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[sS][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[sS]", "[BCDGHKMNRS@V@Y]");
            //System.out.println("here?V " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[vV][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[vV]", "[ABCDGHKMNRSV@W@Y]");
            //System.out.println("here?W " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[wW][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[wW]", "[ABDHKMNRTUVW@Y]");
            //System.out.println("here?Y " + regexp.length());
            //if (regexp.matches("[a-z[A-Z]]*(?<!@)[yY][a-z[A-Z]]*"))
                regexp = regexp.replaceAll("(?<!@)[yY]", "[BCDHKMNSTUVWY]");
        }
        else {
            regexp = regexp.replaceAll("(?<!@)[bB]", "[CGTU]");
            regexp = regexp.replaceAll("(?<!@)[dD]", "[AGTU]");
            regexp = regexp.replaceAll("(?<!@)[hH]", "[ACTU]");
            regexp = regexp.replaceAll("(?<!@)[kK]", "[GTU]");
            regexp = regexp.replaceAll("(?<!@)[mM]", "[AC]");
            regexp = regexp.replaceAll("(?<!@)[nN]", "[ACGTU]");
            regexp = regexp.replaceAll("(?<!@)[rR]", "[AG]");
            regexp = regexp.replaceAll("(?<!@)[sS]", "[CG]");
            regexp = regexp.replaceAll("(?<!@)[vV]", "[ACG]");
            regexp = regexp.replaceAll("(?<!@)[wW]", "[ATU]");
            regexp = regexp.replaceAll("(?<!@)[yY]", "[CTU]");
            degen = true;
        }

        //System.out.println("here?# " + regexp.length());
        if (regexp.indexOf("#") != -1)
            regexp = regexp.replaceAll("#", "[ABCDGHKMNRSTUVWY]*?");

        //Unprotect
        regexp = regexp.replaceAll("@", "");

        //System.out.println("This Heap Space?");

        return regexp;
    }


     /**
     * Determines the GC content of a string of nucleotides, returns a double
     * in the range from 0 to 1
     *
     */
    public double[] gcContent() {
        String seq = sequence;
        seq = seq.toUpperCase();
        int len = seqLength();
        double gcMin = 0;
        double gcMax = 0;
        double[] gc = new double[2];
        for (int i = 0; i < len; i++) {
            String n = seq.substring(i, i + 1);
            if (n.matches("[CGS]")) {
                gcMin++;
                gcMax++;
            }
            else if (n.matches("[RYMKBDHVN]")) {
                gcMax++;
            }
        }
        gcMin = gcMin / len;
        gcMax = gcMax / len;

        gc[0] = gcMin;
        gc[1] = gcMax;
        return gc;
    }


    /**
     * Determines the approximate melting point (Celsius) of a sequence of DNA
     * using a Nearest-Neighbor method, assuming 1.0 M [NaCl] and
     * 50 nM [primer].
     *
     */
    public double meltingTemp () {
        
        /* Resources:
         * http://en.wikipedia.org/wiki/DNA_melting#Nearest-neighbor_method
         * http://www.basic.northwestern.edu/biotools/oligocalc.html
         * http://dna.bio.puc.cl/cardex/servers/dnaMATE/tm-pred.html
         */

        String seq = sequence;
        int len = seqLength();
        double concP = 50 * java.lang.Math.pow(10, -9);
        double dH = 0;
        double dS = 0;
        double logCt = 0;
        double R = 1.987;
        double temp;
        String pair;
        seq = seq.toUpperCase();

        // Checks terminal base pairs
        char init = seq.charAt(0);
        if (init == 'G' || init == 'C') {
            dH += 0.1;
            dS += -2.8;
        }
        else if (init == 'A' || init == 'T') {
            dH += 2.3;
            dS += 4.1;
        }
        init = seq.charAt(len - 1);
        if (init == 'G' || init == 'C') {
            dH += 0.1;
            dS += -2.8;
        }
        else if (init == 'A' || init == 'T') {
            dH += 2.3;
            dS += 4.1;
        }

        // Checks nearest neighbor pairs
        for (int i = 0; i < len - 1; i++) {
            pair = seq.substring(i,i+2);
            if (pair.equals("AA") || pair.equals("TT")) {
                dH += -7.9;
                dS += -22.2;
            }
            else if (pair.equals("AG") || pair.equals("CT")) {
                dH += -7.8;
                dS += -21.0;
            }
            else if (pair.equals("AT")) {
                dH += -7.2;
                dS += -20.4;
            }
            else if (pair.equals("AC") || pair.equals("GT") ) {
                dH += -8.4;
                dS += -22.4;
            }
            else if (pair.equals("GA") || pair.equals("TC")) {
                dH += -8.2;
                dS += -22.2;
            }
            else if (pair.equals("GG") || pair.equals("CC")) {
                dH += -8.0;
                dS += -19.9;
            }
            else if (pair.equals("GC")) {
                dH += -9.8;
                dS += -24.4;
            }
            else if (pair.equals("TA")) {
                dH += -7.2;
                dS += -21.3;
            }
            else if (pair.equals("TG") || pair.equals("CA")) {
                dH += -8.5;
                dS += -22.7;
            }
            else if (pair.equals("CG") ) {
                dH += -10.6;
                dS += -27.2;
            }
        }

        // Checks for self-complementarity
        int mid;
        if (len % 2 == 0) {
            mid = len / 2;
            if (seq.substring(0, mid).equals(new NucSeq(seq.substring(mid,len)).revComp())) {
                dS += -1.4;
            }
        }
        else {
            mid = (len - 1) / 2;
            if (seq.substring(0, mid).equals(new NucSeq(seq.substring(mid + 1,len)).revComp())) {
                dS += -1.4;
            }
        }

        // dH is in kCal, dS is in Cal - equilibrating units
        dH = dH * 1000;

        // logCt = java.lang.Math.log(1 / concP);
        logCt = java.lang.Math.log(concP);

        temp = (dH / (dS + (R * logCt))) - 273.15;

        //return temp;
        return temp;
    }


    /**
     * Reverse complement this NucSeq.  The case will be saved with this
     * operation, and the annotations will be repositioned.
     */
    public void revCompThis() {
        char currentchar;
        StringBuffer seq = new StringBuffer();
        boolean[] newLower = new boolean[ sequence.length() ];

        for ( int x = (sequence.length() - 1); x >= 0; x-- ) {
            currentchar = sequence.charAt( x );
            char appendChar = ' ';
            switch ( currentchar ) {   // (Assume N is an integer variable.)
                case 'A':
                    if ( isRNA ) {
                        appendChar = 'U';
                    } else {
                        appendChar = 'T';
                    }
                    break;
                case 'T':
                    appendChar = 'A';
                    break;
                case 'C':
                    appendChar = 'G';
                    break;
                case 'G':
                    appendChar = 'C';
                    break;
                case '&':
                    appendChar = '&';
                    break;
                case 'R':
                    appendChar = 'Y';
                    break;
                case 'Y':
                    appendChar = 'R';
                    break;
                case 'M':
                    appendChar = 'K';
                    break;
                case 'K':
                    appendChar = 'M';
                    break;
                case 'W':
                    appendChar = 'W';
                    break;
                case 'S':
                    appendChar = 'S';
                    break;
                case 'B':
                    appendChar = 'V';
                    break;
                case 'D':
                    appendChar = 'H';
                    break;
                case 'H':
                    appendChar = 'D';
                    break;
                case 'V':
                    appendChar = 'B';
                    break;
                case 'N':
                    appendChar = 'N';
                    break;
                case 'U':
                    appendChar = 'A';
                    break;
                default:
                    break;
            }
            seq.append( appendChar );
            if ( lowerArray[x] ) {
                appendChar = Character.toLowerCase( appendChar );
                newLower[sequence.length() - x - 1] = true;
            } else {
                newLower[sequence.length() - x - 1] = false;
            }
        }

        //update the sequence
        lowerArray = newLower;
        if(changeSeq(seq.toString())) {
            //Invert all annotations
            for ( Annotation an : annotations ) {
                an.invert( seq.length() );
            }
        }

    }

    /**
     * Subroutine revComp returns the reverse complement of theSequence,
     * in all uppercase as a String.  To actually reverse complement this
     * NucSeq, and also invert its annotations, use revCompThis instead.
     *
     * @return a String that is the reverse complement
     */
    public String revComp() {
        StringBuffer seq = new StringBuffer();
        char currentchar;

        for ( int x = (sequence.length() - 1); x >= 0; x-- ) {
            currentchar = sequence.charAt( x );
            char outchar = ' ';
            switch ( currentchar ) {   // (Assume N is an integer variable.)
                case 'A':
                    if ( isRNA ) {
                        outchar = 'U';
                    } else {
                        outchar = 'T';
                    }
                    break;
                case 'T':
                    outchar = 'A';
                    break;
                case 'C':
                    outchar = 'G';
                    break;
                case 'G':
                    outchar = 'C';
                    break;
                case '&':
                    outchar = '&';
                    break;
                case 'R':
                    outchar = 'Y';
                    break;
                case 'Y':
                    outchar = 'R';
                    break;
                case 'M':
                    outchar = 'K';
                    break;
                case 'K':
                    outchar = 'M';
                    break;
                case 'W':
                    outchar = 'W';
                    break;
                case 'S':
                    outchar = 'S';
                    break;
                case 'B':
                    outchar = 'V';
                    break;
                case 'D':
                    outchar = 'H';
                    break;
                case 'H':
                    outchar = 'D';
                    break;
                case 'V':
                    outchar = 'B';
                    break;
                case 'N':
                    outchar = 'N';
                    break;
                case 'U':
                    outchar = 'A';
                    break;
                default:
                    break;
            }
            seq.append( outchar );
        }  // end for loop
        return seq.toString();
    }

    /**
     * Get the NucSeq in Genbank format from the user currently logged in
     * @return a String in Genbank format
     */
    public String getGenbank() {
        //return getGenbank( new Person[]{ Collector.getCurrentUser() } );
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Get the NucSeq in Genbank format with annotations from a specific
     * list of users
     * @param users
     * @return a String in Genbank format
     */
    public String getGenbank( Person[] users ) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Annotate the NucSeq from the complete database worth of features
     * If Person is null, uses the current user from Collector
     *
     * The Person is used as the author of the Annotation
     */
    public void autoAnnotate( Person user) {
        System.out.println( "I'm autoannotating your nucSeq from all database features" );
        autoAnnotate( null, user, false );
    }

    /**
     * Annotate the NucSeq from the features in a particular Collection
     * If Person is null, uses the current user from Collector
     *
     * The Person is used as the author of the Annotation
     * @param col
     */
    public void autoAnnotate( Collection col, Person user ) {
        Set<Feature> allfeatures = new HashSet<Feature>(col.getAll(Feature.class));
        System.out.println( "Autoannotating with all features from a particular collection:" );
        for ( Feature f : allfeatures ) {
            System.out.println( "autoannotate with " + f.getId() );
        }
        autoAnnotate( allfeatures, user, true );
    }

    /**
     * Relay for other two autoAnnotate calls, but can also use this in Apps directly.
     * For using all features in the database, set constrainTo == false.
     *
     * To search a particularly Collection,
     * call autoAnnotate ( Collection col, Person user ).  To search all collections of
     * a particular user call autoAnnotate ( Person user )
     *
     * The Person is used as the author of the Annotation
     *
     * @param onlyFeatures the list of Feature UUIDs requested for autoannotation
     * @param user  the Person to be set as author of the Annotation
     * @param constrainTo true if should constrain annotations to supplied list, otherwise false
     */
    public void autoAnnotate( Set<Feature> onlyFeatures, Person user, Boolean constrainTo ) {
        if ( !featuresInitiated ) {
            initiateFeatureTable();
        }

        String revcomp = this.revComp();

        for ( Feature f : featureTable.keySet()) {
            if(constrainTo) {
                if ( onlyFeatures != null ) {
                    if ( !onlyFeatures.contains( f ) ) {
                        System.out.println( f.getId() + " is not requested" );
                        continue;
                    }
                }
            }

            try {
                testFeature(featureTable.get(f), f, user, revcomp );
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Iterated method that compares one feature sequence to the target
     *
     * @param teststring
     * @param uuid
     * @param user
     * @param revcomp
     */
    private void testFeature(String teststring, Feature f, Person user, String revcomp ) {
        Pattern p = Pattern.compile(teststring);

        //Check Feature exact matches in forward orientation
        String[] text = {sequence};
        for (int i = 0; i < text.length; i++) {
            Matcher matcher = p.matcher(text[i]);
            while (matcher.find()) {
                System.out.println( "start=" + matcher.start() + " end = " + matcher.end());
                if(f==null || f.isDeleted()) {
                    return;
                }
                int start = matcher.start();
                int end = matcher.end();
                if(f.isCDS()) {
                    try {
                        //For CDS features, if the 5' sequences is a start codon, include that in annotation
                        String fiveprime = text[i].substring(start-3, start);
                        System.out.println("fiveprime is " + fiveprime);
                        if(fiveprime.equals("ATG") || fiveprime.equals("TTG") || fiveprime.equals("GTG")) {
                            start = start-3;
                        }
                    } catch(Exception e) {
                    }
                    try {
                        //For CDS features, if the 3' sequences is a stop codon, include that in annotation
                        String threeprime = text[i].substring(end, end+3);
                        System.out.println("threeprime is " + threeprime);
                        if(threeprime.equals("TAA") || threeprime.equals("TGA") || threeprime.equals("TAG")) {
                            end = end+3;
                        }
                    } catch(Exception e) {
                    }
                }
                Annotation annot = new Annotation( f, this, null, null, start, end, user, true, null );
                System.out.println( "I found a forward annotation at " + start );
                //setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.ANNOTATION_TO_NUCSEQ);
            }

        }

        //Check it as reverse complement
        String[] text2 = {revcomp};
        for (int i = 0; i < text.length; i++) {
            Matcher matcher = p.matcher(text2[i]);
            while (matcher.find()) {
                System.out.println( "start=" + matcher.start() + " end = " + matcher.end());
                if(f==null || f.isDeleted()) {
                    return;
                }
                int index = sequence.length() - matcher.start();
                int start = index - teststring.length();
                int end = index;
                if(f.isCDS()) {
                    try {
                        //For CDS features, if the 5' sequences is a an RC stop codon, include it
                        String fiveprime = text[i].substring(start-3, start);
                        System.out.println("fiveprime is " + fiveprime);
                        if(fiveprime.equals("TTA") || fiveprime.equals("TCA") || fiveprime.equals("CTA")) {
                            start = start-3;
                        }
                    } catch(Exception e) {
                    }
                    try {
                        //For CDS features, if the 3' sequences is a an RC start codon, include that in annotation
                        String threeprime = text[i].substring(end, end+3);
                        System.out.println("threeprime is " + threeprime);
                        if(threeprime.equals("CAT") || threeprime.equals("CAA") || threeprime.equals("CAC")) {
                            end = end+3;
                        }
                    } catch(Exception e) {
                    }
                }
                Annotation annot = new Annotation( f, this, null, null, start, end, user, false, null );
                System.out.println( "I found a reverse annotation at " + start );
                //setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.ANNOTATION_TO_NUCSEQ);
            }
        }
    }

    public void removeAnnotations() {
        annotations = new HashSet<Annotation>();
        //setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.ANNOTATION_TO_NUCSEQ);
    }

    void setLocked( boolean isit ) {
        isLocked = isit;
    }

    /**
     * Add a user-defined non-Feature Annotation
     * @param annotation
     */
    void addAnnotation(Annotation annotation){
        annotations.add(annotation);
    }

    public ArrayList<Integer> find( NucSeq seq ) {
        ArrayList<Integer> out = new ArrayList<Integer>();
        String testSeq = seq.sequence;

        int start = 0;
        searchforward:
        while ( true ) {
            String test = sequence.substring( start );
            System.out.println( test.substring( 0, 100 ) );
            int b = test.indexOf( testSeq );
            if ( b > 0 ) {
                out.add( b + start );
                System.out.println( "for adding: " + b );
                start += b;
                start += 1;
                System.out.println( "New start: " + start );
            } else {
                break searchforward;
            }
        }

        testSeq = seq.revComp();
        start = 0;
        searchreverse:
        while ( true ) {
            String test = sequence.substring( start );
            System.out.println( test.substring( 0, 100 ) );
            int b = test.indexOf( testSeq );
            if ( b > 0 ) {
                out.add( b + start );
                System.out.println( "for adding: " + b );
                start += b;
                start += 1;
                System.out.println( "New start: " + start );
            } else {
                break searchreverse;
            }
        }
        if ( out.size() == 0 ) {
            System.out.println( "I didn't find any" );
        }
        return out;
    }

    public String translate( int frame ) {
        return translate( frame, sequence.length() );
    }

    public String translate( int start, int end ) {
        int extra = (end - start) % 3;
        if ( extra % 3 > 0 ) {
            System.out.println( "You gave me an invalid translation query: " + sequence.substring( start, end ) );
            return "*";
        }
        String seq = sequence.substring( start, end );

        int value;
        int i = 0;
        String acodon = "";
        String outSeq = "";
        while ( i < seq.length() ) {
            acodon = seq.substring( i, i + 3 );
            if ( translation.containsKey( acodon ) ) {
                int anum = translation.get( acodon );
                outSeq += (char) anum;
            } else {
                outSeq += "?";
            }
            i = i + 3;
        }

        return outSeq;
    }

    public char getCharAt( int i ) {
        return sequence.charAt( i );
    }

    public int seqLength() {
        return sequence.length();
    }

    public boolean isDegenerate() {
        return isDegenerate;
    }

    public boolean isRNA() {
        return isRNA;
    }

    public boolean isLocked() {
        return isLocked;
    }

    public boolean isCircular() {
        return isCircular;
    }

    public boolean isSingleStranded() {
        return isSingleStranded;
    }


    public Set<ObjectId> getAnnotationLinks() {
        Set<ObjectId> out = new HashSet<ObjectId>();
        for (Annotation a : annotations){
            out.add(a.getId());
        }
        return out;
    }

    /**
     * Returns the user-formatted version of the String
     * @return
     */
    @Override
    public String toString() {
        StringBuffer seq = new StringBuffer();
        for ( int i = 0; i < sequence.length(); i++ ) {
            char letter = sequence.charAt( i );
            if ( lowerArray[i] ) {
                letter = Character.toLowerCase( letter );
            }
            seq.append( letter );
        }
        return seq.toString();
    }

    /**
     * Returns the unformatted all-caps string, used
     * for bioinformaticcy treatment
     */
    public String getSeq() {
        return sequence;
    }

    /**
     * Returns the unformatted all-caps string with all
     * degeneracy positions replaced by regex
     * @return
     */
    public String getMatcher() {
        String out = sequence;
        out.replaceAll("N", ".");
        out.replaceAll("R", ".");
        out.replaceAll("K", ".");
        out.replaceAll("S", ".");
        return out;
    }

    
    //TODO: proper composition in Part, Vector, feature, oligo
    /**
     * Change the sequence of this NucSeq.  Parts, Vectors,
     * features, and oligos "lock" their NucSeq...you must
     * call changeSeq from the part, vector, Feature, or oligo
     * to change their sequence.
     *
     * @param newseq
     */
    public boolean changeSeq( String newseq ) {
        if ( isLocked ) {
            return false;
        }
        //ADD UNDO HERE FOR THESEQUENCE AND ANNOTATIONS, THEN CLEAR ANNOTATIONS
        return APIchangeSeq( newseq );
    }

    boolean APIchangeSeq( String newseq ) {
        return initiateNucSeq(newseq);
    }

    private static final ImmutableMap<String, Integer> translation = ImmutableMap.<String, Integer>builder()
            .put( "TTT",  70 ) 
       .put( "TTC",  70  )
        .put( "TTA",  76  )
        .put( "TTG",  76  )
        .put( "CTT",  76  )
        .put( "CTC",  76  )
        .put( "CTA",  76  )
        .put( "CTG",  76  )
        .put( "ATT",  73  )
        .put( "ATC",  73  )
        .put( "ATA",  73  )
        .put( "ATG",  77  )
        .put( "GTT",  86  )
        .put( "GTC",  86  )
        .put( "GTA",  86  )
        .put( "GTG",  86  )
        .put( "TCT",  83  )
        .put( "TCC",  83  )
        .put( "TCA",  83  )
        .put( "TCG",  83  )
        .put( "CCT",  80  )
        .put( "CCC",  80  )
        .put( "CCA",  80  )
        .put( "CCG",  80  )
        .put( "ACT",  84  )
        .put( "ACC",  84  )
        .put( "ACA",  84  )
        .put( "ACG",  84  )
        .put( "GCT",  65  )
        .put( "GCC",  65  )
        .put( "GCA",  65  )
        .put( "GCG",  65  )
        .put( "TAT",  89  )
        .put( "TAC",  89  )
        .put( "TAA",  42  )
        .put( "TAG",  42  )
        .put( "CAT",  72  )
        .put( "CAC",  72  )
        .put( "CAA",  81  )
        .put( "CAG",  81  )
        .put( "AAT",  78  )
        .put( "AAC",  78  )
        .put( "AAA",  75  )
        .put( "AAG",  75  )
        .put( "GAT",  68  )
        .put( "GAC",  68  )
        .put( "GAA",  69  )
        .put( "GAG",  69  )
        .put( "TGT",  67  )
        .put( "TGC",  67  )
        .put( "TGA",  42  )
        .put( "TGG",  87  )
        .put( "CGT",  82  )
        .put( "CGC",  82  )
        .put( "CGA",  82  )
        .put( "CGG",  82  )
        .put( "AGT",  83  )
        .put( "AGC",  83  )
        .put( "AGA",  82  )
        .put( "AGG",  82  )
        .put( "GGT",  71  )
        .put( "GGC",  71  )
        .put( "GGA",  71  )
        .put( "GGG",  71  )
            .build();

    public static void refreshFeatureTable() {
        generateFeatureTable( false, true );
    }

    public static void initiateFeatureTable() {
        generateFeatureTable( true, false );
    }

    private static void generateFeatureTable( boolean init, boolean backgroundMode ) {
        //FIXME
       /* Feature[] features = Persistor.get().get(Feature.class);

        for ( int i = 0; i < features.length; i++ ) {
            try {
            String astring = features[i].getSequence().getSeq();
            String string2 = astring.toUpperCase();
            featureTable.put(features[i], string2.replaceAll("N", "."));
            } catch(Exception e) {
                featureTable.put(features[i], "XXXXXXXXXXXXXXXXXX");
            }
        }
        featuresInitiated = true;*/
    }

    /**
     * Relayed from feature constructors to add the local memory feature to the table of autoannotations
     * @param afeature
     */
    static void addFeatureToTable(Feature afeature) {
        if(afeature==null) {
            return;
        }
        String seq = afeature.getSequence().getSeq();
        featureTable.put(afeature, seq.replaceAll("N", "."));
    }

    /**
     * This is the general method called to perform all biofafety checks
     * on a DNA sequence.  Called from part, vector, and Feature factory
     * methods.
     *
     * @return the biosafety level of this NucSeq
     */
    /*Short performBiosafetyCheck() {
       System.out.println( "performBiosafetyCheck triggered" );
       short rg = -1;
       rg = getBSLfromServer();
       //If it's RG3+, show a special message
       if ( rg == 4 ) {
           ImageIcon bslicon = ImageUtilities.loadImageIcon( "org/clothocore/images/BIOHAZARD.png", false );
           JOptionPane.showMessageDialog( null,
                                          "You have executed a risk group check on a sequence that came back\n"
                   + "as Risk Group 4.\nSuch materials could be highly dangerous!\n"
                   + "You should examine your design closer before proceeding.",
                                          "Risk Group 4 material!",
                                          JOptionPane.INFORMATION_MESSAGE,
                                          bslicon );
       }
       if ( rg == 5 ) {
           ImageIcon bslicon = ImageUtilities.loadImageIcon( "org/clothocore/images/BIOHAZARD.png", false );
           JOptionPane.showMessageDialog( null,
                                          "You have executed a risk group check on a sequence that came back\n"
                   + "as being highly similar to a select agent.\nSuch materials could be highly dangerous and potential illegal!\n"
                   + "You should examine your design closer before proceeding.",
                                          "Select Agent detected!",
                                          JOptionPane.INFORMATION_MESSAGE,
                                          bslicon );
       }
       return rg;
   }

   /**
    * This queries the bsl server and parses the xml to
    * get the biosafety level of the NucSeq
    *
    * @return the biosafety level of this NucSeq
    */
   /*private short getBSLfromServer(){
        //If it's already failed 3 times, then don't bother anymore
        if(failCount>3) {
            return -1;
        }

        //Form the URL query
        String seq = getSeq();
        URL urlRobot;
	try {
            String urlstr = _BSLServerURL + "\"" + seq + "\"";
	    urlRobot = new URL(urlstr);
	} catch (Exception e) {
	    e.printStackTrace();
            return -1;
	}

        XMLParser myParser = null;
	try {
            //Starts reading the URL
	    InputStream urlRobotStream = urlRobot.openStream();

            //Read the file and access it in an xmlParser, then close it
            try {
                myParser = new XMLParser(urlRobotStream, "output" );
            } catch(Exception e) {
                System.out.println("Biosafety server data could not be parsed");
                failCount++;
                updateBSLServer(_updateBSLURL1);
                return -1;
            }
	    urlRobotStream.close();
	} catch (java.net.ConnectException e) {
	    System.out.println("Biosafety server not available");
            //Keeps count of failures.  Once pass 3 give up.
            failCount++;
            updateBSLServer(_updateBSLURL1);
            return -1;
	} catch (java.io.IOException ex) {
            System.out.println("Biosafety server not available");
            failCount++;
            updateBSLServer(_updateBSLURL1);
            return -1;
        }

        if(myParser==null) {
            System.out.println("Biosafety server information could not be parsed");
            failCount++;
            updateBSLServer(_updateBSLURL1);
            return-1;
        }

        try {
            String bslvalue = myParser.getFirstTag("bsl");
            short bsl = Short.parseShort(bslvalue);
            System.out.println("Biosafety server returning risk group " + bsl);
            if(bsl == (short) 0) {
                bsl=1;
            }
            return bsl;
        } catch(Exception e) {
            failCount++;
            updateBSLServer(_updateBSLURL1);
           return -1;
        }
   }

   /**
    * The first time this class is called, set the biosafety BLAST
    * server from preferences.  If no preference is set, try retrieving
    * it from an XML file online
    */
   /*static {
       String tempurl = Collator.getPreference("NucSeqBSLServerAddress");
       try {
           URL url = new URL(tempurl);
       } catch (Exception ex) {
           updateBSLServer("http://www.bu.edu/ece-clotho/xmlfes/updatebsl.xml");
       }
   }

   /**
    * Request that NucSeq update its biosafety server.
    */
    /*public static void updateBSLServer(String url) {
	//Form URL of the file
	URL urlRobot = null;
	try {
	    urlRobot = new URL(url);
	} catch (Exception e) {
	    e.printStackTrace();
            if(!url.equals(_updateBSLURL2)) {
                updateBSLServer(_updateBSLURL2);
                return;
            }
	}

        XMLParser myParser = null;
	try {
            //Starts reading the URL
	    InputStream urlRobotStream = urlRobot.openStream();

            //Read the file and access it in an xmlParser, then close it
            myParser = new XMLParser(urlRobotStream, "update" );
	    urlRobotStream.close();
	} catch (Exception e) {
	    e.printStackTrace();
            if(!url.equals(_updateBSLURL2)) {
                updateBSLServer(_updateBSLURL2);
                return;
            }
	}
        if(myParser==null) {
            return;
        }

        String newurl = myParser.getFirstTag("url");
        try {
            URL urly = new URL(newurl);
        } catch (MalformedURLException ex) {
            return;
        }

        Collator.putPreference("NucSeqBSLServerAddress", newurl);
        _BSLServerURL = newurl;
        System.out.println("The new BSL server from " + urlRobot.getPath() + " is " + newurl);
    }

/*-----------------
    variables
-----------------*/
    /*private static final String _updateBSLURL1 = "http://www.bu.edu/ece-clotho/xmlfiles/updatebsl.xml";
    private static final String _updateBSLURL2 = "http://andersonlab.qb3.berkeley.edu/Software/updatebsl.xml";
    private static short failCount = 0;
    private static String _BSLServerURL = "http://cidar1.bu.edu/cgi-bin/tst.pl?";*/

    protected boolean[] lowerArray;

   /* public static class NucSeqDatum extends ObjBaseDatum {

        public boolean isDegenerate = false;   // if it has N's or R's and so on
        public boolean isRNA = false;          // if it has U's
        public boolean isLinear = false;          // if it '.''s
        public boolean isSingleStranded = false;   // if its an oligo
        public boolean isCircular = false;     // if its a plasmid
        public String theSequence;
        public HashSet<String> annotations = new HashSet<String>();  //The list of annoations
        public HashSet<String> removeAnnotations = new HashSet<String>(); // ? What is is this?
        public boolean isLocked = false;

    }*/
    
    //Is there a better way to cache all features?
    //should we cache all features?
    private static boolean featuresInitiated = false;
    private static boolean initiating = false;
    
    
    //converts feature sequence from N to . when loaded
    //uppercases feature sequence
   private static Map<Feature, String> featureTable;

}
