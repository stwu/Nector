package Bio;

import java.io.*;
import java.util.*;
public class DnaSequenceFasta extends SequenceFasta {
	private static char [] complement = { 'T', 'B', 'G', 'D', 'E', 'F', 'C',
	                               'H', 'I', 'J', 'K', 'L', 'M', 'N',
							       'O', 'P', 'Q',      'R', 'S', 'A',
								   'A', 'V', 'W',      'X', 'Y', 'Z'};
    private static char[] capital = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J','K', 'L', 'M', 'N',
                       'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};

    public DnaSequenceFasta(){
        super();
	}
    public DnaSequenceFasta(String description, String sequence){
        super(description, sequence);
	}
    public DnaSequenceFasta(DnaSequenceFasta fasta){
        super(fasta); 
	}
	public boolean sequenceCheck(){
	   char [] seq = super.getSequence().toCharArray();
	   int [] count = new int[26];
	   int i, len;
	   for(i = 0; i < 26; i++){
          count[i] = 0;
	   }
	   len = seq.length;
       for(char a: seq){
		  count[a - 'A']++; 
	   }
       int sum = 0;
	   char [] tmp = {'A', 'C', 'G', 'T', 'N'};
	   for(char t: tmp){
		   sum += count[t - 'A'];   
	   }
	   if(sum == len){
	      return true;
	   }
	   else{
		  System.err.println("Error! Sequences should only contains A, C, G, T, N");
          return false;
	   }
	}
	public boolean descriptionCheck(){
       if(!super.getDescription().substring(0,1).equals(">")){
		    System.err.println("Error! Description should begin with \">\" in FASTA format!");
            return false;
	   }
	   else{
            return true;
	   }
	}
	public boolean formatCheck(){
	   //check quality area
       if(!descriptionCheck()){
            return false;
	   }
	   //check sequence
       if(!sequenceCheck()){
            return false;
	   }
	   return true;
	}
	public DnaSequenceFasta complementaryReverse(){
        int i, j, len, mid;
		int left, right;
		char [] array;
		len = super.getLength();
		array = super.getSequence().toCharArray();
		//System.out.println(super.getSequence());
		//System.out.println(len);
		mid = (len / 2) + (len % 2);
		for(i = 0, j = len - 1; i < mid; i++, j--){
           left = array[i];
		   right = array[j];
           array[i] = complement[right - 'A']; 
           array[j] = complement[left - 'A']; 
		}
		return new DnaSequenceFasta(super.getDescription(), new String(array));
	}
	//for testing purpose
	public static void main(String [] args){
           DnaSequenceFasta fasta = new DnaSequenceFasta(">seq1", "aaaaCCCCGGGGGTTTTTNa");
           System.out.print(fasta.formatCheck());
		   DnaSequenceFasta fasta2 = new DnaSequenceFasta(">seq2", "GGGCCCAATTN");
           System.out.print(fasta2);
           System.out.print(fasta2.complementaryReverse());
	}
}
