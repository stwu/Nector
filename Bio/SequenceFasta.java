package Bio;
import java.io.*;
import java.util.*;

public class SequenceFasta implements Comparable<SequenceFasta> {
	private String description;//Desciption of sequences, treat it like ID, don't include ">"
    private String sequence;//Sequence like ACGT

	public SequenceFasta(){
       description = null;
	   sequence = null;
	}
    public SequenceFasta(SequenceFasta fasta){
       setFasta(fasta);
	}

	public SequenceFasta(String description, String sequence){
	   this.description = description;
       this.sequence = sequence.toUpperCase();
	}

	public void setSequence(String sequence){
       this.sequence = sequence.toUpperCase();
	}

	public void setDescription(String description){
       this.description = description;
	}
	public void setFasta(SequenceFasta fasta){
	   description = fasta.getDescription(); 
       sequence = fasta.getSequence();
	}

	public String getSequence(){
       return sequence;
	}

	public String getDescription(){
       return description;
	}
	public int getLength(){
       return sequence.length(); 
	}
	/*@Override
	public SequenceFasta clone() throws CloneNotSupportedException {
       return (SequenceFasta) super.clone();
	}*/
	@Override
	public int compareTo(SequenceFasta other){
       return sequence.length() - other.sequence.length();
	}
	@Override
	public String toString(){
       return description + "\n" + sequence + "\n";
	}
	//The following is for testing purpose
    public static void main(String [] args){
		  SequenceFasta fasta = new SequenceFasta(">seq1", "aaaaCCCCGGGGGTTTTT");
          System.out.print(fasta);
		  SequenceFasta fasta2 = new SequenceFasta(">seq2", "AACCCCGGGGGTTTTT");
          System.out.print(fasta2);
		  SequenceFasta fasta3 = fasta;
          System.out.print(fasta3);
		  System.out.println(fasta.compareTo(fasta2));
		  System.out.println(fasta2.compareTo(fasta));
		  System.out.println(fasta.compareTo(fasta3));
		  System.out.println(fasta.getLength());
		  System.out.println(fasta2.getLength());
		  System.out.println(fasta3.getLength());
	}
}
