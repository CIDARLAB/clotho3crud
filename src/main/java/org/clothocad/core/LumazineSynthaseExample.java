package org.clothocad.core;

import org.clothocad.core.ClothoBuilder;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.nosecurity.NoSecurityModule;
import org.clothocad.core.util.JSON;
import org.clothocad.model.Annotation;
import org.clothocad.model.BasicModule;
import org.clothocad.model.BioDesign;
import org.clothocad.model.CompositeModule;
import org.clothocad.model.Feature;
import org.clothocad.model.Variable;
import org.clothocad.model.Feature.FeatureRole;
import org.clothocad.model.Influence;
import org.clothocad.model.Influence.InfluenceType;
import org.clothocad.model.Module;
import org.clothocad.model.Module.ModuleRole;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;
import org.clothocad.model.Strain;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/*
* This example is based off of an entry from the iGEM registry
* http://parts.igem.org/Part:BBa_K542008
*/
public class LumazineSynthaseExample {
	
    public static void main(String[] args) {
        ClothoBuilder builder = new ClothoBuilder(new NoSecurityModule(), new JongoModule());
        Persistor persistor = builder.get(Persistor.class);

        ObjectId designID = createDesignK542008(persistor);
        
        BioDesign design = (BioDesign) readObjBase(designID, persistor);
        
        updateObjBaseName(design.getModule().getId(), "Nicholas Roehner", persistor);
    }
    
    public static void updateObjBaseName(ObjectId objectID, String rename, Persistor persistor) {
    	ObjBase obj = persistor.get(ObjBase.class, objectID);
    	obj.setName(rename);
    	persistor.save(obj, true);
    }
    
    public static Object readObjBase(ObjectId objectID, Persistor persistor) {
    	return persistor.get(ObjBase.class, objectID);
    }
	
	public static ObjectId createDesignK542008(Persistor persistor) {
		/*
		* Authors
		*/
		
		Person nic = new Person("Nic");
		Person anonymous = new Person("Anonymous");
		Person anthony = new Person("Anthony Vuong");
		Person vinay = new Person("Vinay S. Mahajan");
		Person roxanne = new Person("Roxanne Shank");
		Person randy = new Person("Randy Rettberg");
		Person reshma = new Person("Reshma Shetty");
		Person jkm = new Person("J. K. M.");
		Person june = new Person("June Rhee");
		Person lisza = new Person("Lisza Bruder");
		Person kristen = new Person("Kristen DeCelle");
		Person adam = new Person("Adam Smith");
		
		/*
		* Sequences
		*/
		
		// Elowitz RBS sequence
		Sequence seqB0034 = new Sequence("B0034 Sequence", "aaagaggagaaa", vinay);

		// Terminator sequences
		Sequence seqB0010 = new Sequence("B0010 Sequence", "ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctc", randy);
		Sequence seqB0012 = new Sequence("B0012 Sequence", "tcacactggctcaccttcgggtgggcctttctgcgtttata", reshma);
		
		// Lumazine Synthase gene sequences
		Sequence seqR0010 = new Sequence("R0010 Sequence", "caatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaag"
				+ "cgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcgga"
				+ "taacaatttcacaca", anonymous);
		Sequence seqK249002 = new Sequence("K249002 Sequence", "atgcagatttatgaaggcaaactgaccgcggaaggcctgcgctttggcattgtggcgagccgctttaaccatgcgc"
				+ "tggtggatcgcctggtggaaggcgcgattgattgcattgtgcgccatggtggtcgcgaagaagatattaccctggtgcgcgtgccgggcagctgggaaattccggtgg"
				+ "cggcgggcgaactggcgcgcaaagaagatattgatgcggtgattgcgattggcgtgctgattgaaggcgcggaaccgcattttgattatattgcgagcgaagtgagca"
				+ "aaggcctggcgaacctgagcctggaactgcgcaaaccgattacctttggcgtgattaccgcggatgaactggaagaagcgattgaacgcgcgggcaccaaacatggca"
				+ "acaaaggctgggaagcggcgctgagcgcgattgaaatggcgaacctgtttaaaagcctgcgctag", roxanne);
		
		// TetR gene sequences
		Sequence seqI13453 = new Sequence("I13453 Sequence", "acattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctgac"
				+ "gctttttatcgcaactctctactgtttctccataccgtttttttgggctagc", jkm);
		Sequence seqC0040 = new Sequence("C0040 Sequence", "atgtccagattagataaaagtaaagtgattaacagcgcattagagctgcttaatgaggtcggaatcgaaggtttaacaa"
				+ "cccgtaaactcgcccagaagctaggtgtagagcagcctacattgtattggcatgtaaaaaataagcgggctttgctcgacgccttagccattgagatgttagataggca"
				+ "ccatactcacttttgccctttagaaggggaaagctggcaagattttttacgtaataacgctaaaagttttagatgtgctttactaagtcatcgcgatggagcaaaagt"
				+ "acatttaggtacacggcctacagaaaaacagtatgaaactctcgaaaatcaattagcctttttatgccaacaaggtttttcactagagaatgcattatatgcactcag"
				+ "cgctgtggggcattttactttaggttgcgtattggaagatcaagagcatcaagtcgctaaagaagaaagggaaacacctactactgatagtatgccgccattattacg"
				+ "acaagctatcgaattatttgatcaccaaggtgcagagccagccttcttattcggccttgaattgatcatatgcggattagaaaaacaacttaaatgtgaaagtgggtc"
				+ "cgctgcaaacgacgaaaactacgctttagtagcttaataacactgatagtgctagtgtagatcac", june);
		
		// ECFP/EYFP gene sequences
		Sequence seqR0040 = new Sequence("R0040 Sequence", "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac", june);
		Sequence seqK331002 = new Sequence("K331002 Sequence", "atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtgaacggcca"
				+ "caagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgt"
				+ "gaccaccctgacctggggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcac"
				+ "catcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggac"
				+ "ggcaacatcctggggcacaagctggagtacaactacatcagccacaacgtctatatcaccgccgacaagcagaagaacggcatcaaggccaacttcaagatccgccaca"
				+ "acatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccg"
				+ "ccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagcactagagccg"
				+ "ccgccgccgccgccgccgccgccgctaa", lisza);
		Sequence seqK249006 = new Sequence("K249006 Sequence", "atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggcc"
				+ "acaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgt"
				+ "gaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcacc"
				+ "atcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacg"
				+ "gcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaac"
				+ "atcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccc"
				+ "tgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagcactagagccgccgc"
				+ "cgccgccgccgccgccgccgctaa", roxanne);
		
		// Double terminator sequences
		// seqB0015 contains seqB0010 and seqB0012
		Sequence seqB0015 = new Sequence("B0015 Sequence", "ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctct"
				+ "ctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata", reshma);
		
		// Lumazine Synthase gene sequences
		// seqJ04500 contains seqR0010 and seqB0034
		Sequence seqJ04500 = new Sequence("J04500 Sequence", "caatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaa"
				+ "gcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggat"
				+ "aacaatttcacacatactagagaaagaggagaaa", kristen);
		// seqK542000 contains seqK249002 and seqB0015
		Sequence seqK542000 = new Sequence("K542000 Sequence", "atgcagatttatgaaggcaaactgaccgcggaaggcctgcgctttggcattgtggcgagccgctttaaccatgcgc"
				+ "tggtggatcgcctggtggaaggcgcgattgattgcattgtgcgccatggtggtcgcgaagaagatattaccctggtgcgcgtgccgggcagctgggaaattccggtg"
				+ "gcggcgggcgaactggcgcgcaaagaagatattgatgcggtgattgcgattggcgtgctgattgaaggcgcggaaccgcattttgattatattgcgagcgaagtgagc"
				+ "aaaggcctggcgaacctgagcctggaactgcgcaaaccgattacctttggcgtgattaccgcggatgaactggaagaagcgattgaacgcgcgggcaccaaacatggc"
				+ "aacaaaggctgggaagcggcgctgagcgcgattgaaatggcgaacctgtttaaaagcctgcgctagtactagagccaggcatcaaataaaacgaaaggctcagtcga"
				+ "aagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata", anthony);
		
		// TetR gene sequence
		// seqP0440 contains seqB0034 and seqC0040 and seqB0015
		Sequence seqP0440 = new Sequence("P0440 Sequence", "aaagaggagaaatactagatgtccagattagataaaagtaaagtgattaacagcgcattagagctgcttaatgaggtc"
				+ "ggaatcgaaggtttaacaacccgtaaactcgcccagaagctaggtgtagagcagcctacattgtattggcatgtaaaaaataagcgggctttgctcgacgccttagc"
				+ "cattgagatgttagataggcaccatactcacttttgccctttagaaggggaaagctggcaagattttttacgtaataacgctaaaagttttagatgtgctttactaag"
				+ "tcatcgcgatggagcaaaagtacatttaggtacacggcctacagaaaaacagtatgaaactctcgaaaatcaattagcctttttatgccaacaaggtttttcactaga"
				+ "gaatgcattatatgcactcagcgctgtggggcattttactttaggttgcgtattggaagatcaagagcatcaagtcgctaaagaagaaagggaaacacctactactga"
				+ "tagtatgccgccattattacgacaagctatcgaattatttgatcaccaaggtgcagagccagccttcttattcggccttgaattgatcatatgcggattagaaaaaca"
				+ "acttaaatgtgaaagtgggtccgctgcaaacgacgaaaactacgctttagtagcttaataacactgatagtgctagtgtagatcactactagagccaggcatcaaata"
				+ "aaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctg"
				+ "cgtttata", randy);
		
		// ECFP/EYFP gene sequences
		// seqK331025 contains seqB0034 and seqK331002
		Sequence seqK331025 = new Sequence("K331025 Sequence", "aaagaggagaaatactagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctg"
				+ "gacggcgacgtgaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgc"
				+ "ccgtgccctggcccaccctcgtgaccaccctgacctggggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccg"
				+ "aaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctg"
				+ "aagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacatcagccacaacgtctatatcaccgccgacaagcagaagaacggcat"
				+ "caaggccaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccga"
				+ "caaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcat"
				+ "ggacgagctgtacaagcactagagccgccgccgccgccgccgccgccgccgctaa", adam);
		// seqK331023 contains seqB0034 and seqK249006
		Sequence seqK331023 = new Sequence("K331023 Sequence", "aaagaggagaaatactagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctg"
				+ "gacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcc"
				+ "cgtgccctggcccaccctcgtgaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccga"
				+ "aggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaa"
				+ "gggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatca"
				+ "aggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgac"
				+ "aaccactacctgagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcat"
				+ "ggacgagctgtacaagcactagagccgccgccgccgccgccgccgccgccgctaa", adam);
			
		// ECFP/EYFP gene sequences
		// seqK331033 contains seqR0040 and seqK331025
		Sequence seqK331033 = new Sequence("K331033 Sequence", "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcactactagagaaagaggagaaata"
				+ "ctagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtgaacggccacaagttcagcgtgtccggcgagggcg"
				+ "agggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctggggcgtgcagtg"
				+ "cttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaacta"
				+ "caagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctgga"
				+ "gtacaactacatcagccacaacgtctatatcaccgccgacaagcagaagaacggcatcaaggccaacttcaagatccgccacaacatcgaggacggcagcgtgcagct"
				+ "cgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaa"
				+ "gcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagcactagagccgccgccgccgccgccgccgccgccg"
				+ "ctaa", adam);
		// seqK331035 contains seqK331023 and seqB0015
		Sequence seqK331035 = new Sequence("K331035 Sequence", "aaagaggagaaatactagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctgga"
				+ "cggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgt"
				+ "gccctggcccaccctcgtgaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccgaagg"
				+ "ctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaaggg"
				+ "catcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaag"
				+ "gtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaac"
				+ "cactacctgagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacg"
				+ "agctgtacaagcactagagccgccgccgccgccgccgccgccgccgctaatactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgtt"
				+ "ttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata", adam);
			
		// Lumazine Synthase gene sequence, TetR gene sequence
		// seqK542001 contains seqJ04500 and seqK542000
		Sequence seqK542001 = new Sequence("K542001 Sequence", "caatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaa"
				+ "gcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggat"
				+ "aacaatttcacacatactagagaaagaggagaaatactagatgcagatttatgaaggcaaactgaccgcggaaggcctgcgctttggcattgtggcgagccgctttaac"
				+ "catgcgctggtggatcgcctggtggaaggcgcgattgattgcattgtgcgccatggtggtcgcgaagaagatattaccctggtgcgcgtgccgggcagctgggaaattc"
				+ "cggtggcggcgggcgaactggcgcgcaaagaagatattgatgcggtgattgcgattggcgtgctgattgaaggcgcggaaccgcattttgattatattgcgagcgaagt"
				+ "gagcaaaggcctggcgaacctgagcctggaactgcgcaaaccgattacctttggcgtgattaccgcggatgaactggaagaagcgattgaacgcgcgggcaccaaacat"
				+ "ggcaacaaaggctgggaagcggcgctgagcgcgattgaaatggcgaacctgtttaaaagcctgcgctagtactagagccaggcatcaaataaaacgaaaggctcagtc"
				+ "gaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata", anthony);
		// seqK542003 contains seqI13453 and seqP0440
		Sequence seqK542003 = new Sequence("K542003 Sequence", "acattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctga"
				+ "cgctttttatcgcaactctctactgtttctccataccgtttttttgggctagctactagagaaagaggagaaatactagatgtccagattagataaaagtaaagtgat"
				+ "taacagcgcattagagctgcttaatgaggtcggaatcgaaggtttaacaacccgtaaactcgcccagaagctaggtgtagagcagcctacattgtattggcatgtaaaa"
				+ "aataagcgggctttgctcgacgccttagccattgagatgttagataggcaccatactcacttttgccctttagaaggggaaagctggcaagattttttacgtaataac"
				+ "gctaaaagttttagatgtgctttactaagtcatcgcgatggagcaaaagtacatttaggtacacggcctacagaaaaacagtatgaaactctcgaaaatcaattagcc"
				+ "tttttatgccaacaaggtttttcactagagaatgcattatatgcactcagcgctgtggggcattttactttaggttgcgtattggaagatcaagagcatcaagtcgct"
				+ "aaagaagaaagggaaacacctactactgatagtatgccgccattattacgacaagctatcgaattatttgatcaccaaggtgcagagccagccttcttattcggcctt"
				+ "gaattgatcatatgcggattagaaaaacaacttaaatgtgaaagtgggtccgctgcaaacgacgaaaactacgctttagtagcttaataacactgatagtgctagtgt"
				+ "agatcactactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacac"
				+ "tggctcaccttcgggtgggcctttctgcgtttata", anthony);
				
		// Lumazine Synthase gene and TetR gene sequence, ECFP/EYFP gene sequence
		// seqK542004 contains seqK542001 and seqK542003
		Sequence seqK542004 = new Sequence("K542004 Sequence", "caatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactgga"
				+ "aagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagc"
				+ "ggataacaatttcacacatactagagaaagaggagaaatactagatgcagatttatgaaggcaaactgaccgcggaaggcctgcgctttggcattgtggcgagccgc"
				+ "tttaaccatgcgctggtggatcgcctggtggaaggcgcgattgattgcattgtgcgccatggtggtcgcgaagaagatattaccctggtgcgcgtgccgggcagctg"
				+ "ggaaattccggtggcggcgggcgaactggcgcgcaaagaagatattgatgcggtgattgcgattggcgtgctgattgaaggcgcggaaccgcattttgattatattg"
				+ "cgagcgaagtgagcaaaggcctggcgaacctgagcctggaactgcgcaaaccgattacctttggcgtgattaccgcggatgaactggaagaagcgattgaacgcgcg"
				+ "ggcaccaaacatggcaacaaaggctgggaagcggcgctgagcgcgattgaaatggcgaacctgtttaaaagcctgcgctagtactagagccaggcatcaaataaaac"
				+ "gaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtt"
				+ "tatatactagagacattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctgacgctttttatcgcaactc"
				+ "tctactgtttctccataccgtttttttgggctagctactagagaaagaggagaaatactagatgtccagattagataaaagtaaagtgattaacagcgcattagagct"
				+ "gcttaatgaggtcggaatcgaaggtttaacaacccgtaaactcgcccagaagctaggtgtagagcagcctacattgtattggcatgtaaaaaataagcgggctttgct"
				+ "cgacgccttagccattgagatgttagataggcaccatactcacttttgccctttagaaggggaaagctggcaagattttttacgtaataacgctaaaagttttagat"
				+ "gtgctttactaagtcatcgcgatggagcaaaagtacatttaggtacacggcctacagaaaaacagtatgaaactctcgaaaatcaattagcctttttatgccaacaa"
				+ "ggtttttcactagagaatgcattatatgcactcagcgctgtggggcattttactttaggttgcgtattggaagatcaagagcatcaagtcgctaaagaagaaagggaa"
				+ "acacctactactgatagtatgccgccattattacgacaagctatcgaattatttgatcaccaaggtgcagagccagccttcttattcggccttgaattgatcatatg"
				+ "cggattagaaaaacaacttaaatgtgaaagtgggtccgctgcaaacgacgaaaactacgctttagtagcttaataacactgatagtgctagtgtagatcactactaga"
				+ "gccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttc"
				+ "gggtgggcctttctgcgtttata", anthony);
		// seqK542005 contains seqK331033 and seqK331035
		Sequence seqK542005 = new Sequence("K542005 Sequence", "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcactactagagaaagaggagaaata"
				+ "ctagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtgaacggccacaagttcagcgtgtccggcgagggcga"
				+ "gggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctggggcgtgcagtg"
				+ "cttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaact"
				+ "acaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctgg"
				+ "agtacaactacatcagccacaacgtctatatcaccgccgacaagcagaagaacggcatcaaggccaacttcaagatccgccacaacatcgaggacggcagcgtgcag"
				+ "ctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgag"
				+ "aagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagcactagagccgccgccgccgccgccgccgccg"
				+ "ccgctaatactagagaaagaggagaaatactagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggc"
				+ "cacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccct"
				+ "cgtgaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagc"
				+ "gcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaag"
				+ "gaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagat"
				+ "ccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagct"
				+ "accagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaag"
				+ "cactagagccgccgccgccgccgccgccgccgccgctaatactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgtt"
				+ "gtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata", anthony);

		// Lumazine Synthase gene and TetR gene and ECFP/EYFP gene sequence
		// seqK542008 contains seqK542004 and seqK542005
		Sequence seqK542008 = new Sequence("Sequence K542008", "caatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaa"
				+ "agcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcgg"
				+ "ataacaatttcacacatactagagaaagaggagaaatactagatgcagatttatgaaggcaaactgaccgcggaaggcctgcgctttggcattgtggcgagccgcttt"
				+ "aaccatgcgctggtggatcgcctggtggaaggcgcgattgattgcattgtgcgccatggtggtcgcgaagaagatattaccctggtgcgcgtgccgggcagctgggaa"
				+ "attccggtggcggcgggcgaactggcgcgcaaagaagatattgatgcggtgattgcgattggcgtgctgattgaaggcgcggaaccgcattttgattatattgcgagc"
				+ "gaagtgagcaaaggcctggcgaacctgagcctggaactgcgcaaaccgattacctttggcgtgattaccgcggatgaactggaagaagcgattgaacgcgcgggcacc"
				+ "aaacatggcaacaaaggctgggaagcggcgctgagcgcgattgaaatggcgaacctgtttaaaagcctgcgctagtactagagccaggcatcaaataaaacgaaaggc"
				+ "tcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttatatact"
				+ "agagacattgattatttgcacggcgtcacactttgctatgccatagcatttttatccataagattagcggatcctacctgacgctttttatcgcaactctctactgt"
				+ "ttctccataccgtttttttgggctagctactagagaaagaggagaaatactagatgtccagattagataaaagtaaagtgattaacagcgcattagagctgcttaat"
				+ "gaggtcggaatcgaaggtttaacaacccgtaaactcgcccagaagctaggtgtagagcagcctacattgtattggcatgtaaaaaataagcgggctttgctcgacgc"
				+ "cttagccattgagatgttagataggcaccatactcacttttgccctttagaaggggaaagctggcaagattttttacgtaataacgctaaaagttttagatgtgctt"
				+ "tactaagtcatcgcgatggagcaaaagtacatttaggtacacggcctacagaaaaacagtatgaaactctcgaaaatcaattagcctttttatgccaacaaggttttt"
				+ "cactagagaatgcattatatgcactcagcgctgtggggcattttactttaggttgcgtattggaagatcaagagcatcaagtcgctaaagaagaaagggaaacacct"
				+ "actactgatagtatgccgccattattacgacaagctatcgaattatttgatcaccaaggtgcagagccagccttcttattcggccttgaattgatcatatgcggatt"
				+ "agaaaaacaacttaaatgtgaaagtgggtccgctgcaaacgacgaaaactacgctttagtagcttaataacactgatagtgctagtgtagatcactactagagccag"
				+ "gcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtg"
				+ "ggcctttctgcgtttatatactagagtccctatcagtgatagagattgacatccctatcagtgatagagatactgagcactactagagaaagaggagaaatactaga"
				+ "tggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtgaacggccacaagttcagcgtgtccggcgagggcgagggc"
				+ "gatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctggggcgtgcagtgctt"
				+ "cagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactaca"
				+ "agacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggag"
				+ "tacaactacatcagccacaacgtctatatcaccgccgacaagcagaagaacggcatcaaggccaacttcaagatccgccacaacatcgaggacggcagcgtgcagct"
				+ "cgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgaga"
				+ "agcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagcactagagccgccgccgccgccgccgccgccgc"
				+ "cgctaatactagagaaagaggagaaatactagatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggc"
				+ "cacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccct"
				+ "cgtgaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagc"
				+ "gcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaag"
				+ "gaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaaga"
				+ "tccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctga"
				+ "gctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtac"
				+ "aagcactagagccgccgccgccgccgccgccgccgccgctaatactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatct"
				+ "gttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata", anthony);
		
		/*
		* Parts
		*/
		
		// Elowitz RBS part
		Part partB0034 = new Part("BBa_B0034", seqB0034, vinay);
		
		// Terminator parts
		Part partB0010 = new Part("BBa_B0010", seqB0010, randy);
		Part partB0012 = new Part("BBa_B0012", seqB0012, reshma);
		
		// Lumazine Synthase gene parts
		Part partR0010 = new Part("BBa_R0010", seqR0010, anonymous);
		Part partK249002 = new Part("BBa_K249002", seqK249002, roxanne);
		
		// TetR gene parts
		Part partI13453 = new Part("BBa_I13453", seqI13453, jkm);
		Part partC0040 = new Part("BBa_C0040", seqC0040, june);
		
		// ECFP/EYFP gene parts
		Part partR0040 = new Part("BBa_R0040", seqR0040, june);
		Part partK331002 = new Part("BBa_K331002", seqK331002, lisza);
		Part partK249006 = new Part("BBa_K249006", seqK249006, roxanne);
		
		// Double terminator part
		// partB0015 is assembled from partB0010 and partB0012
		Part partB0015 = new Part("BBa_B0015", seqB0015, reshma);
		partB0015.createAssembly().setParts(Arrays.asList(partB0010, partB0012));
		
		// Lumazine Synthase gene parts
		// partJ04500 is assembled from partR0010 and partB0034
		Part partJ04500 = new Part("BBa_J04500", seqJ04500, kristen);
		partJ04500.createAssembly().setParts(Arrays.asList(partR0010, partB0034));
		// partK542000 is assembled from partK249002 and partB0015
		Part partK542000 = new Part("BBa_K542000", seqK542000, anthony);
		partK542000.createAssembly().setParts(Arrays.asList(partK249002, partB0015));
		
		// TetR gene part
		Part partP0440 = new Part("BBa_P0440", seqP0440, randy);
		
		// ECFP/EYFP gene parts
		// partK331025 is assembled from partB0034 and partK331002
		Part partK331025 = new Part("BBa_K331025", seqK331025, adam);
		partK331025.createAssembly().setParts(Arrays.asList(partB0034, partK331002));
		// partK331023 is assembled from partB0034 and partK249006
		Part partK331023 = new Part("BBa_K331023", seqK331023, adam);
		partK331023.createAssembly().setParts(Arrays.asList(partB0034, partK249006));
		
		// ECFP/EYFP gene parts
		// partK331033 is assembled from partR0040 and partK331025
		Part partK331033 = new Part("BBa_K331033", seqK331033, adam);
		partK331033.createAssembly().setParts(Arrays.asList(partR0040, partK331025));
		// partK331035 is assembled from partK331023 and partB0015
		Part partK331035 = new Part("BBa_K331035", seqK331035, adam);
		partK331035.createAssembly().setParts(Arrays.asList(partK331023, partB0015));
		
		// Lumazine Synthase gene part, TetR gene part
		// partK542001 is assembled from partJ04500 and partK542000
		Part partK542001 = new Part("BBa_K542001", seqK542001, anthony);
		partK542001.createAssembly().setParts(Arrays.asList(partJ04500, partK542000));
		// partK542003 is assembled from partI13453 and partP0440
		Part partK542003 = new Part("BBa_K542003", seqK542003, anthony);
		partK542003.createAssembly().setParts(Arrays.asList(partI13453, partP0440));
		
		// Lumazine Synthase gene and TetR gene part, ECFP/EYFP gene part
		// partK542004 is assembled from partK542001 and partK542003
		Part partK542004 = new Part("BBa_K542004", seqK542004, anthony);
		partK542004.createAssembly().setParts(Arrays.asList(partK542001, partK542003));
		// partK542005 is assembled from partK331033 and partK331035
		Part partK542005 = new Part("BBa_K542005", seqK542005, anthony);
		partK542005.createAssembly().setParts(Arrays.asList(partK331033, partK331035));
		
		// Lumazine Synthase gene and TetR gene and ECFP/EYFP gene part
		// partK542008 is assembled from partK542004 and partK542005
		Part partK542008 = new Part("BBa_K542008", seqK542008, anthony);
		partK542008.createAssembly().setParts(Arrays.asList(partK542004, partK542005));
		
		/*
		* Features and Annotations
		*/
		
		// Elowitz RBS feature
		Feature featB0034 = constructFeature("RBS (Elowitz 1999) -- defines RBS efficiency", FeatureRole.RBS, seqB0034, vinay);
		seqB0034.createAnnotation("Conserved", 5, 8, true, vinay);
		
		// Terminator features
		Feature featB0010 = constructFeature("T1 from E. coli rrnB", FeatureRole.TERMINATOR, seqB0010, randy);
		seqB0010.createAnnotation("stem loop", 12, 55, true, randy);
		Feature featB0012 = constructFeature("TE from coliphageT7", FeatureRole.TERMINATOR, seqB0012, reshma);
		seqB0012.createAnnotation("stem loop", 8, 27, true, reshma);
		seqB0012.createAnnotation("Stop", 34, 34, true, reshma);
		seqB0012.createAnnotation("PolyA", 28, 41, true, reshma);
		
		// Lumazine Synthase gene features
		Feature featR0010 = constructFeature("promoter (lacI regulated)", FeatureRole.PROMOTER, seqR0010, anonymous);
		constructFeature("end of LacI coding region (inactive)", FeatureRole.CDS, seqR0010, 1, 88, anonymous);
		constructFeature("CAP binding site", FeatureRole.OPERATOR, seqR0010, 89, 126, anonymous);
		constructFeature("-35", FeatureRole.PROMOTER, seqR0010, 137, 142, anonymous);
		constructFeature("-10", FeatureRole.PROMOTER, seqR0010, 161, 166, anonymous);
		constructFeature("LacI binding site", FeatureRole.OPERATOR, seqR0010, 166, 200, anonymous);
		seqR0010.createAnnotation("start", 173, 173, true, anonymous);
		Feature featK249002 = constructFeature("Lumazine Synthase", FeatureRole.CDS, seqK249002, roxanne);
		
		// TetR gene features
		Feature featI13453 = constructFeature("Pbad promoter", FeatureRole.PROMOTER, seqI13453,  jkm);
		Feature featC0040 = constructFeature("tetracycline repressor from transposon Tn10 (+LVA)", FeatureRole.CDS, seqC0040, june);
		constructFeature("tetR", FeatureRole.CDS, seqC0040, 1, 620, june);
		constructFeature("SsrA", FeatureRole.DEGRADATION_TAG, seqC0040, 621, 654, june);
		constructFeature("barcode", FeatureRole.BARCODE, seqC0040, 661, 685, june);
		
		// ECFP/EYFP gene features
		Feature featR0040 = constructFeature("TetR repressible promoter", FeatureRole.PROMOTER, seqR0040, june);
		constructFeature("TetR 1", FeatureRole.OPERATOR, seqR0040, 1, 19, june);
		constructFeature("-35", FeatureRole.PROMOTER, seqR0040, 20, 25, june);
		constructFeature("TetR 2", FeatureRole.OPERATOR, seqR0040, 26, 44, june);
		constructFeature("-10", FeatureRole.PROMOTER, seqR0040, 43, 48, june);
		Feature featK331002 = constructFeature("ECFP with C-terminal Arginine Tag", FeatureRole.CDS, seqK331002, lisza);
		seqK331002.createAnnotation("start", 1, 3, true, lisza);
		constructFeature("Arginine Tag", FeatureRole.LOCALIZATION_TAG, seqK331002, 725, 753, lisza);
		seqK331002.createAnnotation("stop", 754, 756, true, lisza);
		Feature featK249006 = constructFeature("Fusion EYFP with C-terminal Arginine Tag", FeatureRole.CDS, seqK249006, roxanne);
		
		// Double terminator feature and annotations
		// seqB0015 is annotated with featB0010 and featB0012
		Feature featB0015 = constructFeature("double terminator (B0010-B0012)", FeatureRole.TERMINATOR, seqB0015, reshma);
		annotateWithFeature(seqB0015, featB0010, 1, 80, reshma);
		annotateWithFeature(seqB0015, featB0012, 89, 129, reshma);
		
		// Lumazine Synthase gene annotations
		// seqJ04500 is annotated with featR0010 and featB0034
		annotateWithFeature(seqJ04500, featR0010, 1, 200, kristen);
		annotateWithFeature(seqJ04500, featB0034, 209, 220, kristen);
		// seqK542000 is annotated with featK249002 and featB0015
		annotateWithFeature(seqK542000, featK249002, 1, 465, anthony);
		annotateWithFeature(seqK542000, featB0015, 474, 602, anthony);
		
		// TetR gene annotations
		// seqP0440 is annotated with featB0034 and featC0040 and featB0015
		annotateWithFeature(seqP0440, featB0034, 1, 12, randy);
		annotateWithFeature(seqP0440, featC0040, 19, 678, randy);
		annotateWithFeature(seqP0440, featB0015, 712, 840, randy);
		
		// ECFP/EYFP gene annotations
		// seqK331025 is annotated with featB0034 and featK331002
		annotateWithFeature(seqK331025, featB0034, 1, 12, adam);
		annotateWithFeature(seqK331025, featK331002, 19, 774, adam);
		// seqK331023 is annotated with featB0034 and featK249006
		annotateWithFeature(seqK331023, featB0034, 1, 12, adam);
		annotateWithFeature(seqK331023, featK249006, 19, 774, adam);
		
		// ECFP/EYFP gene annotations
		// seqK331033 is annotated with featR0040 and featB0034 and featK331002
		annotateWithFeature(seqK331033, featR0040, 1, 54, adam);
		annotateWithFeature(seqK331033, featB0034, 63, 74, adam);
		annotateWithFeature(seqK331033, featK331002, 81, 836, adam);
		// seqK331035 is annotated with featB0034 and featK249006 and featB0015
		annotateWithFeature(seqK331035, featB0034, 1, 12, adam);
		annotateWithFeature(seqK331035, featK249006, 19, 774, adam);
		annotateWithFeature(seqK331035, featB0015, 783, 911, adam);
		
		// Lumazine Synthase gene feature and annotations, TetR gene feature and annotations
		Feature featK542001 = constructFeature("pLacI Regulated Lumazine Synthase", FeatureRole.GENE, seqK542001, anthony);
		// seqK542001 is annotated with featR0010 and featB0034 and featK249002 and featB0015
		annotateWithFeature(seqK542001, featR0010, 1, 200, anthony);
		annotateWithFeature(seqK542001, featB0034, 209, 220, anthony);
		annotateWithFeature(seqK542001, featK249002, 227, 691, anthony);
		annotateWithFeature(seqK542001, featB0015, 700, 828, anthony);
		Feature featK542003 = constructFeature("pBAD Regulated TetR", FeatureRole.GENE, seqK542003, anthony);
		// seqK542003 is annotated with featI13453 and featB0034 and featC0040 and featB0015
		annotateWithFeature(seqK542003, featI13453, 1, 130, anthony);
		annotateWithFeature(seqK542003, featB0034, 139, 150, anthony);
		annotateWithFeature(seqK542003, featC0040, 157, 816, anthony);
		annotateWithFeature(seqK542003, featB0015, 850, 978, anthony);
		
		// Lumazine Synthase gene and TetR gene annotations, ECFP/EYFP gene feature and annotations
		// seqK542004 is annotated with featK542001 and featK542003
		annotateWithFeature(seqK542004, featK542001, 1, 828, anthony);
		annotateWithFeature(seqK542004, featK542003, 837, 1814, anthony);
		Feature featK542005 = constructFeature("pTet Regulated Arg-tagged ECFP and EYFP (FRET Reporter)", FeatureRole.GENE, seqK542005, anthony);
		// seqK542005 is annotated with featR0040 and featB0034 and featK331002 and featK249006 and featB0015
		annotateWithFeature(seqK542005, featR0040, 1, 54, anthony);
		annotateWithFeature(seqK542005, featB0034, 63, 74, anthony);
		annotateWithFeature(seqK542005, featK331002, 81, 836, anthony);
		annotateWithFeature(seqK542005, featB0034, 845, 856, anthony);
		annotateWithFeature(seqK542005, featK249006, 863, 1618, anthony);
		annotateWithFeature(seqK542005, featB0015, 1627, 1755, anthony);
		
		// Lumazine Synthase gene and TetR gene and ECFP/EYFP gene annotations
		// seqK542008 is annotated with featK542001 and featK542003 and featK542005
		annotateWithFeature(seqK542008, featK542001, 1, 828, anthony);
		annotateWithFeature(seqK542008, featK542003, 837, 1814, anthony);
		annotateWithFeature(seqK542008, featK542005, 1823, 3577, anthony);
	
		/*
		* Modules
		*/
		
		// RBS module
		Module modB0034 = constructBasicModule("Translation via Elowitz RBS", ModuleRole.TRANSLATION, featB0034, nic);
		
		// Terminator modules
		Module modB0010 = constructBasicModule("Transcription with T1 from E. coli rrnB", ModuleRole.TRANSCRIPTION, featB0010, nic);
		Module modB0012 = constructBasicModule("Transcription with TE from coliphageT7", ModuleRole.TRANSCRIPTION, featB0012, nic);
		
		// Lumazine Synthase gene modules
		Module modR0010 = constructBasicModule("Transcription via pLac", ModuleRole.TRANSCRIPTION, featR0010, nic);
		Module modK249002 = constructBasicModule("Lumazine Synthase Expression", ModuleRole.EXPRESSION, featK249002, nic);
		
		// TetR gene modules
		Module modI13453 = constructBasicModule("Transcription via pBad", ModuleRole.TRANSCRIPTION, featI13453, nic);
		Module modC0040 = constructBasicModule("TetR Expression", ModuleRole.EXPRESSION, featC0040, nic);
		
		// ECFP/EYFP gene modules
		Module modR0040 = constructBasicModule("Transcription via pTet", ModuleRole.TRANSCRIPTION, featR0040, nic);
		Module modK331002 = constructBasicModule("ECFP Expression", ModuleRole.EXPRESSION, featK331002, nic);
		Module modK249006 = constructBasicModule("EYFP Expression", ModuleRole.EXPRESSION, featK249006, nic);
		
		// Double terminator module
		Module modB0015 = constructCompositeModule("Transcription with Double Terminator", ModuleRole.TRANSCRIPTION, modB0010, modB0012, nic);
		
		// Lumazine Synthase gene modules
		// modJ04500 is comprised of modR0010 and modB0034
		Module modJ04500 = constructCompositeModule("Expression via pLac and Elowitz RBS", ModuleRole.EXPRESSION, modR0010, modB0034, nic);
		// modK542000 is comprised of modK249002 and modB0015
		Module modK542000 = constructCompositeModule("Lumazine Synthase Expression", ModuleRole.EXPRESSION, modK249002, modB0015, nic);
		
		// TetR gene module
		// modP0440 is comprised of modB0034 and modC0040 and modB0015
		CompositeModule modP0440 = constructCompositeModule("TetR Expression", ModuleRole.EXPRESSION, modB0034, modC0040, nic);
		modP0440.addSubModule(modB0015);
		
		// ECFP/EYFP gene modules
		// modK331025 is comprised of modB0034 and modK331002
		Module modK331025 = constructCompositeModule("ECFP Expression", ModuleRole.EXPRESSION, modB0034, modK331002, nic);
		// modK331023 is comprised of modB0034 and modK249006
		Module modK331023 = constructCompositeModule("EYFP Expression", ModuleRole.EXPRESSION, modB0034, modK249006, nic);
		
		// ECFP/EYFP gene modules
		// modK331033 is comprised of modR0040 and modK331025
		Module modK331033 = constructCompositeModule("ECFP Expression", ModuleRole.EXPRESSION, modR0040, modK331025, nic);
		// modK331035 is comprised of modK331023 and modB0015
		Module modK331035 = constructCompositeModule("EYFP Expression", ModuleRole.EXPRESSION, modK331023, modB0015, nic);
		
		// Lumazine Synthase gene module, TetR gene module
		// modK542001 is comprised of modJ04500 and modK542000
		Module modK542001 = constructCompositeModule("Lumazine Synthase Expression", ModuleRole.EXPRESSION, modJ04500, modK542000, nic);
		// modK542003 is comprised of modI13453 and modP0440
		Module modK542003 = constructCompositeModule("TetR Expression", ModuleRole.EXPRESSION, modI13453, modP0440, nic);
		
		// Lumazine Synthase gene and TetR gene module, ECFP/EYFP gene module
		// modK542004 is comprised of modK542001 and modK542003
		Module modK542004 = constructCompositeModule("Lumazine Microcompartment Controller", ModuleRole.COMPARTMENTALIZATION, modK542001, modK542003, 
				nic);
		// modK542005 is comprised of modK331033 and modK331035
		Module modK542005 = constructCompositeModule("FRET Reporter", ModuleRole.REPORTER, modK331033, modK331035, nic);
		
		// Lumazine Synthase gene and TetR gene and ECFP/EYFP gene module and influence
		// modK542008 is comprised of modK542004 and modK542005
		Module modK542008 = constructCompositeModule("FRET Colocalization", ModuleRole.LOCALIZATION, modK542004, modK542005, nic);
		Influence inflC0040_R0040 = new Influence("Repression via TetR", featC0040, featR0040, InfluenceType.REPRESSION, nic);
		Set<Influence> inflsK542008 = new HashSet<Influence>();
		inflsK542008.add(inflC0040_R0040);
		modK542008.setInfluences(inflsK542008);
		
		/*
		* Strain
		*/
		
		Strain strainEcoliDH5a = new Strain("E. coli DH5alpha", nic);
		
		/*
		* Part Designs
		*/
		
		// Elowitz RBS design and parameter
		BioDesign desB0034 = constructDesign("RBS (Elowitz 1999) -- defines RBS efficiency", "RBS based on Elowitz repressilator.", 
				partB0034, modB0034, vinay);
		Variable translEff = new Variable("Translation Efficiency", vinay);
		desB0034.createParameter(1, translEff);
		
		// Terminator designs and parameters
		BioDesign desB0010 = constructDesign("T1 from E. coli rrnB", "Transcriptional terminator consisting of a 64 bp stem-loop.", 
				partB0010, modB0010, randy);
		BioDesign desB0012 = constructDesign("TE from coliphageT7", "Transcription terminator for the E.coli RNA polymerase. (This "
				+ "is a bad terminator, it is a promoter in the reverse direction.)", partB0034, modB0034, reshma);
		Variable forTermEff = new Variable("Forward Termination Efficiency", reshma);
		Variable revTermEff = new Variable("Reverse Termination Efficiency", reshma);
		desB0012.createParameter(0.309, forTermEff);
		desB0012.createParameter(-0.368, revTermEff);
		
		// Lumazine Synthase gene designs
		BioDesign desR0010 = constructDesign("promoter (lacI regulated)", "This part is an inverting regulator sensitive to LacI and CAP. "
				+ "It contains two protein binding sites."
				+ " The first binds the CAP protein, which is generally present in E.coli and is asocciated with cell health and "
				+ "availability of glucose. The second binds LacI protein. In the absence of LacI protein and CAP protein, this part "
				+ "promotes transcription. In the presence of LacI protein and CAP protein, this part inhibits transcription. LacI "
				+ "can be inhibited by IPTG. LacI is coded by BBa_C0010. This is a direct copy of bases 0365739 through 0365540 of "
				+ "the E. coli K-12 MG1655 genome, Genbank NC_000913 in reverse complement form. It is the natural promoter for "
				+ "the LacZYA operon. It includes the tail end of the LacI gene coding region, but no promoter region for that "
				+ "partial gene. ", partR0010, modR0010, anonymous);
		BioDesign desK249002 = constructDesign("Lumazine Synthase", "Lumazine Synthase is an enzyme which creates Lumazine, a product "
				+ "which aggregates forming a "
				+ "hollow spheroid which can act as a mirocompartment, or artificial organelle. The Lumazine forms negatively "
				+ "charged pores, which can be used to introduce proteins. The proteins which are being introduced into the "
				+ "microcompartment must be equipped with an Arginine Tag. Seebeck, F. P., Woycechowsky, K. J., Zhuang, W., "
				+ "Rabe, J. P., and Hilvert, D., (2006). A simple tagging system for protein encapulation. J. Am. Chem. Soc. 128, "
				+ "4516-4517.", partK249002, modK249002, roxanne);
		
		// TetR gene designs
		BioDesign desI13453 = constructDesign("Pbad promoter", "PBad promoter from I0500 without AraC. Has been used as a second promoter "
				+ "in a system containing BBa_I0500 (PBad+AraC). "
				+ "In this system, it showed behavior qualitatively indistinguishable from the BBa_I0500 copy of PBad. Has not been "
				+ "tested independent of AraC. A second part, BBa_I13458, should allow decoupling of PBad and AraC. ", 
				partI13453, modI13453, jkm);
		BioDesign desC0040 = constructDesign("tetracycline repressor from transposon Tn10 (+LVA)", "Coding region for the TetR protein "
				+ "without the Ribosome Binding Site. Modified with an LVA tail "
				+ "for rapid degradation of the protein and faster fall time for the emission. TetR binds to the pTet regulator "
				+ "(Part:BBa_R0040). aTc (anhydrotetracycline) binds to TetR and inhibits its operation.", 
				partC0040, modC0040, june);
		
		// ECFP/EYFP gene designs
		BioDesign desR0040 = constructDesign("TetR repressible promoter", "Sequence for pTet inverting regulator. Promoter is "
				+ "constitutively ON and repressed by TetR. TetR "
				+ "repression is inhibited by the addition of tetracycline or its analog, aTc.", 
				partR0040, modR0040, june);
		BioDesign desK331002 = constructDesign("ECFP with C-terminal Arginine Tag", "This protein is the ECFP fused with the "
				+ "C-terminal Arginine tag found in BBa_K249005. It is "
				+ "designed to be part of a construct for targeting into the lumazine synthase microcompatment (BBa_K331000).", 
				partK331002, modK331002, lisza);
		BioDesign desK249006 = constructDesign("Fusion EYFP with C-terminal Arginine Tag", "This Yellow Fluorescent protein has "
				+ "been fused with 10 C-terminal Arginines to target it "
				+ "into the Lumazine-based microcompartment.", partK249006, modK249006, roxanne);
		
		// Double terminator design and parameters
		// desB0015 is comprised of desB0010 and desB0012
		BioDesign desB0015 = constructDesign("double terminator (B0010-B0012)", "Double terminator consisting of BBa_B0010 and BBa_B0012. This "
				+ "is the most commonly used terminator. It seems to be reliable. Note, however, that Part:BBa_B0014 is a better part for "
				+ "forward and reverse termination.", partB0015, modB0015, reshma);
		desB0015.addSubDesign(desB0010);
		desB0015.addSubDesign(desB0012);
		desB0015.createParameter(0.984, forTermEff);
		desB0015.createParameter(0.295, revTermEff);
		
		/*
		* Intermediate Device Designs
		*/

		// Lumazine Synthase gene designs
		// desJ04500 is comprised of desR0010 and desB0034
		BioDesign desJ04500 = constructDesign("IPTG inducible promoter with RBS", "R0010 and B0034 will be digested and ligated "
				+ "together in the manner described on the \"Registry of Standard Biological Parts\" website.", 
				partJ04500, modJ04500, kristen);
		desJ04500.addSubDesign(desR0010);
		desJ04500.addSubDesign(desB0034);
		// desK542000 is comprised of desK249002 and desB0015
		BioDesign desK542000 = constructDesign("Lumazine Synthase with Transciptional Terminator (No Promoter)", "Made by assembling "
				+ "BBa_K249002 with BBa_B0015. Since this part is lacking the promoter, Lumazine Synthase production "
				+ "may be regulated by the addition different promoters upstream. Regulation of Lumazine Synthase will be "
				+ "dependent on the promoter being "
				+ "utilized. Intermediate part to assemble the \"Co-localization Construct\". (BBa_K542008) This part was "
				+ "characterized in BBa_K542008 in "
				+ "E. coli strain DH5alpha to demonstrate its functionality as an intermediate in that construct.", 
				partK542000, modK542000, strainEcoliDH5a, anthony);
		desK542000.addSubDesign(desK249002);
		desK542000.addSubDesign(desB0015);
		
		// TetR gene designs
		// desP0440 is comprised of desB0034 and desC0040 and desB0015
		BioDesign desP0440 = constructDesign("PoPS -> TetR [S0151]", "-- No description -- ", 
				partP0440, modP0440, randy);
		desP0440.addSubDesign(desB0034);
		desP0440.addSubDesign(desC0040);
		desP0440.addSubDesign(desB0015);
		
		// ECFP/EYFP gene designs
		// desK331025 is comprised of desB0034 and desK331002
		BioDesign desK331025 = constructDesign("RBS with C-terminal Oligo Arginine - ECFP Fusion", "This part has a oligo arginine "
				+ "sequence fused to the C-terminus of an enhanced cyan fluorescent protein. This will be a "
				+ "component of the proof-of-concept part BBa_K331019 which will show localization of proteins into the "
				+ "interior of a microcompartment formed by the assemly of lumazine synthase proteins.", 
				partK331025, modK331025, adam);
		desK331025.addSubDesign(desB0034);
		desK331025.addSubDesign(desK331002);
		// desK331023 is comprised of desB0034 and desK249006
		BioDesign desK331023 = constructDesign("RBS with C-terminal Oligo Arginine - EYFP Fusion", "This part has a oligo arginine "
				+ "sequence fused to the C-terminus of an enhanced yellow fluorescent protein. This will "
				+ "be a component of the proof-of-concept part BBa_K331019 which will show localization of proteins into the "
				+ "interior of a microcompartment formed by the assemly of lumazine synthase proteins.", 
				partK331023, modK331023, adam);
		desK331023.addSubDesign(desB0034);
		desK331023.addSubDesign(desK249006);
		
		// ECFP/EYFP gene designs
		// desK331033 is comprised of desR0040 and desK331025
		BioDesign desK331033 = constructDesign("Tet Repressible C-terminal Arg Tagged ECFP (no terminator)", "This is BBa_K331025 "
				+ "under the control of a tetracycline repressible promoter BBa_R0040. This will be a component of "
				+ "a proof-of-concept part that will show localization of proteins into the interior of a microcompartment "
				+ "formed by the assembly of lumazine synthase proteins.", 
				partK331033, modK331033, adam);
		desK331033.addSubDesign(desR0040);
		desK331033.addSubDesign(desK331025);
		// desK331035 is comprised of desK331023 and desB0015
		BioDesign desK331035 = constructDesign("C-terminal Arg Tagged EYFP with transcriptional terminator", "This is BBa_K331023 "
				+ "with a transcriptional terminator BBa_B0015 added behind. This will be a component of a "
				+ "proof-of-concept part that will show localization of proteins into the interior of a microcompartment "
				+ "formed by the assembly of lumazine synthase proteins.", 
				partK331035, modK331035, adam);
		desK331035.addSubDesign(desK331023);
		desK331035.addSubDesign(desB0015);
		
		/*
		* Device Designs
		*/
		
		// Lumazine Synthase gene design, TetR gene design
		// desK542001 is comprised of desJ04500 and desK542000
		BioDesign desK542001 = constructDesign("pLacI Regulated Lumazine Synthase", "Made by assembling BBa_J04500 with BBa_K542000. "
				+ "IPTG inducible Lumazine Synthase. Intermediate part to assemble the "
				+ "\"Co-localization Construct\". (BBa_K542008) This part was characterized in BBa_K542008 in E. coli strain "
				+ "DH5alpha to demonstrate its functionality as an intermediate in that construct.", 
				partK542001, modK542001, strainEcoliDH5a, anthony);
		desK542001.addSubDesign(desJ04500);
		desK542001.addSubDesign(desK542000);
		// desK542003 is comprised of desI13453 and desP0440
		BioDesign desK542003 = constructDesign("pBAD Regulated TetR", "Made by assembling BBa_I13453 with BBa_P0440. Production of "
				+ "TetR is regulated by the pBAD promoter (BBa_I13453). TetR is "
				+ "an inhibitor of the constitutively \"on\" promoter pTet--BBa_R0040). Therefore, this part may be used in "
				+ "conjunction with the pTet "
				+ "promoter as an \"inverter\"; pTet will be \"turned off\" in the presence of arabinose. (ie. arabinose -> "
				+ "pBAD \"on\" -> TetR -> pTet \"off\") Intermediate part to assemble the \"Co-localization Construct\". "
				+ "(BBa_K542008) This part was characterized in BBa_K542008 in E. coli strain DH5alpha to demonstrate its "
				+ "functionality as an intermediate in that construct.", 
				partK542003, modK542003, strainEcoliDH5a, anthony);
		desK542003.addSubDesign(desI13453);
		desK542003.addSubDesign(desP0440);
		
		// Lumazine Synthase gene and TetR gene design, ECFP/EYFP gene design
		// desK542004 is comprised of desK542001 and desK542003
		BioDesign desK542004 = constructDesign("pLacI Regulated Lumazine Synthase with Transciptional Terminator and pBAD", "Made by "
				+ "assembling BBa_K542001 with BBa_K542003. Note that the pLacI promoter (BBa_R0010) is regulating the Lumazine "
				+ "Synthase gene, not the pBAD promoter (BBa_I13453) which is downstream. Therefore, Lumazine Synthase production "
				+ "is IPTG inducible. The "
				+ "microcompartment portion to the \"Co-localization Construct\". (BBa_K542008) This part was characterized in "
				+ "BBa_K542008 in E. coli strain DH5alpha to demonstrate its functionality as an intermediate in that construct.", 
				partK542004, anthony);
		desK542004.addStrain(strainEcoliDH5a);
		desK542004.addSubDesign(desK542001);
		desK542004.addSubDesign(desK542003);
		// desK542005 is comprised of desK331033 and desK331035
		BioDesign desK542005 = constructDesign("pTet Regulated Arg-tagged ECFP and EYFP (FRET Reporter)", "Made by assembling BBa_K331033 "
				+ "with BBa_K331035. The pTet promoter is constitutively \"on\"; thus, Arg-tagged ECFP and "
				+ "Arg-tagged EYFP will be produced constitutively. "
				+ "By placing this part downstream of BBa_K542003, expression of the fluorescence proteins may be \"turned off\" "
				+ "in the presence of arabinose. Fluorescence/Forster Resonance Energy Transfer (FRET) is a distance-dependent "
				+ "phenomenon in which the excitation of a donor fluorophore leads to emission by an acceptor fluorophore; "
				+ "this occurs if the two fluorophores are within a certain distance to each other. "
				+ "The distance is dependent on the FRET pair used. The emission spectrum of the donor fluorophore must overlap with "
				+ "the excitation spectrum of the acceptor fluorophore for FRET to occur. FRET is explained in further detail on the "
				+ "Lethbridge 2009 Wiki. This phenomenon will allow for characterizing co-localization within the Lumazine Synthase "
				+ "microcompartment. Intermediate part to assemble the \"Co-localization Testing Construct\". (BBa_K542008) This "
				+ "part was characterized in BBa_K542008 in E. coli strain DH5alpha to demonstrate its functionality as an "
				+ "intermediate in that construct.", partK542005, modK542005, strainEcoliDH5a, anthony);
		desK542005.addSubDesign(desK331033);
		desK542005.addSubDesign(desK331035);

		// Lumazine Synthase gene and TetR gene and ECFP/EYFP gene design
		// desK542008 is comprised of desK542004 and desK542005
		BioDesign desK542008 = constructDesign("pLacI Regulated Lumazine Synthase and pBAD Inverse-Regulated Arg-tagged ECFP "
				+ "and EYFP", "Made by assembling BBa_K542004 with BBa_K542005. The intended purpose of the Lumazine "
				+ "Synthase microcompartment device is to specifically target proteins that have been tagged with a positively "
				+ "charged oligopeptide into the cavity formed by the microcompartment. Furthermore, multiple unique proteins "
				+ "may be targeted into the cavity (1) for the purposes of increasing the efficiency of a metabolic pathway, to "
				+ "name just one example. This test construct is designed to demonstrate that proteins with positively charged "
				+ "tags can be localized into the cavity of the microcompartment formed by the oligomerization of Lumazine Synthase "
				+ "monomers. The expression of Lumazine Synthase and the fluorescent proteins are independently regulated by two "
				+ "separate promoters. Lumazine Synthase is regulated by the pLacI promoter and the fluorescent proteins are "
				+ "regulated by the pBAD inverter. Since the two are independently controlled, "
				+ "Lumazine Synthase microcompartments may be formed in the presence or absence of the fluorescent proteins "
				+ "(ie. absence or presence of arabinose, respectively). Because the fluorescent proteins are tagged with a "
				+ "positively-charged poly-arginine tag and the Lumazine Synthase is mutated with a negative interior (BBa_K249002 "
				+ "and Lethbridge 2009 Modeling), enhanced cyan fluorescent protein (ECFP) and "
				+ "enhanced yellow fluorescent protein (EYFP) should be targeted into the microcompartments. ECFP and EYFP are a "
				+ "well known FRET pair. (1) Seebeck, F. P., Woycechowsky, K. J., Zhuang, W., "
				+ "Rabe, J. P., and Hilvert, D., (2006). A simple tagging system for protein encapulation. J. Am. Chem. Soc. 128, "
				+ "4516-4517.", partK542008, modK542008, strainEcoliDH5a, anthony);
		desK542008.addSubDesign(desK542004);
		desK542008.addSubDesign(desK542005);
		
		return persistor.save(desK542008);
	}
	
	/*
	* Convenience construction methods
	*/
	
	public static Feature constructFeature(String name, FeatureRole role, Sequence sequence, Person author) {
		Feature feature = new Feature(name, role, author);
		feature.setSequence(sequence);
		annotateWithFeature(sequence, feature, 1, sequence.getSequence().length(), author);
		return feature;
	}
	
	public static Feature constructFeature(String name, FeatureRole role, Sequence sequence, int start, int end, Person author) {
		Feature feature = new Feature(name, role, author);
		feature.setSequence(new Sequence(name + " Sequence", sequence.getSequence().substring(start - 1, end), author));
		annotateWithFeature(sequence, feature, start, end, author);
		return feature;
	}
	
	public static void annotateWithFeature(Sequence sequence, Feature feature, int start, int end, Person author) {
		Annotation annotation = sequence.createAnnotation(feature.getName(), start, end, true, author);
		annotation.setFeature(feature);
	}
	
	public static BasicModule constructBasicModule(String name, ModuleRole role, Feature feature, Person author) {
		Set<Feature> features = new HashSet<Feature>();
		features.add(feature);
		return new BasicModule(name, role, features, author);
	}
	
	public static CompositeModule constructCompositeModule(String name, ModuleRole role, Module subModule1, Module subModule2, 
			Person author) {
		Set<Module> subModules = new HashSet<Module>();
		subModules.add(subModule1);
		subModules.add(subModule2);
		return new CompositeModule(name, role, subModules, author);
	}
	
	public static BioDesign constructDesign(String name, String description, Part part, Person author) {
		BioDesign design = new BioDesign(name, author);
		Set<Part> parts = new HashSet<Part>();
		parts.add(part);
		design.setParts(parts);
		return design;
	}
	
	public static BioDesign constructDesign(String name, String description, Part part, Module module, Person author) {
		BioDesign design = constructDesign(name, description, part, author);
		design.setModule(module);
		return design;
	}
	
	public static BioDesign constructDesign(String name, String description, Part part, Module module, Strain strain, Person author) {
		BioDesign design = constructDesign(name, description, part, module, author);
		Set<Strain> strains = new HashSet<Strain>();
		strains.add(strain);
		design.setStrains(strains);
		return design;
	}
	
}
