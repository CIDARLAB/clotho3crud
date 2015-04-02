package org.clothocad.model;

import java.io.Serializable;
import java.util.List;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import javax.validation.constraints.NotNull;
import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;

@Data()
@NoArgsConstructor
public class Polynucleotide extends ObjBase implements Serializable {

	//@Getter
	//private ObjectId id;
	//@Getter
	private String description, accession;
	//@Getter
	private boolean isLinear, isSingleStranded;
	//@Getter
	private Date submissionDate;

	//@Getter
	private List<Highlight> highlights;

	@NotNull
	private String sequence;

	/**
	*public Polynucleotide ( String id, String name, String inputSeq,
	*	String description, String accession, String submissionDate,
	*	boolean strandedness, boolean linearity ) {
	*
    *    super(name);
	*
    *    super.setId(new ObjectId(id));
    *    this.sequence = inputSeq;
    *    this.description = description;
    *    this.accession = accession;
	*
    *    //consider doing this with a switch statement too
    *    HashMap<String, Integer> month = new HashMap<String, Integer>();
	*	 month.put("JAN", 1);
	*	 month.put("FEB", 2);
	*	 month.put("MAR", 3);
	*	 month.put("APR", 4);
	*	 month.put("MAY", 5);
	*	 month.put("JUN", 6);
	*	 month.put("JUL", 7);
	*	 month.put("AUG", 8);
	*	 month.put("SEP", 9);
	*	 month.put("OCT", 10);
	*	 month.put("NOV", 11);
	*	 month.put("DEC", 12);
	*
    *    GregorianCalendar date = new GregorianCalendar();
    *    String[] dArr = submissionDate.split("-");
    *    date.set(Integer.parseInt(dArr[2]), month.get(dArr[1]),
    *    	Integer.parseInt(dArr[0]));
    *    this.submissionDate = Date( date.getTime() );
	*
    *    this.isSingleStranded = strandedness;
    *    this.isLinear = linearity;
    *    
    *}
	**/
}
