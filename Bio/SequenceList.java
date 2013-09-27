package Bio;
import java.io.*;
import java.util.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SequenceList{
     private ArrayList<SequenceFasta> sequenceList;
	 private int CAPACITY = 400;
     public SequenceList(){
        sequenceList = new ArrayList<SequenceFasta>();
	 }
	 public SequenceList(String fileName){
		//SequenceList();
        sequenceList = new ArrayList<SequenceFasta>();
		getFromFileFasta(fileName); 
	 }
	 public void getFromFileFasta(String fileName){
		String line, desc = null;
		DnaSequenceFasta fasta = null;
		int i = 0;
		StringBuilder seq = new StringBuilder(CAPACITY);
		seq.setLength(0);
		try{
	      Scanner in = new Scanner(new BufferedReader(new FileReader(fileName)));
	      while(in.hasNextLine()){
		     line = in.nextLine();
			 if(line.substring(0, 1).equals(">")){
			 	if(i > 0){
                    fasta = new DnaSequenceFasta(desc, seq.toString());
					sequenceList.add(fasta);
					seq.setLength(0);
				}
			    desc = new String(line); 	
				i++;
			 }
			 else{
				seq.append(line);
			 }
		  }
		  if(i > 0){
             fasta = new DnaSequenceFasta(desc, seq.toString());
			 sequenceList.add(fasta);
		 	 seq.setLength(0);
		  }
	      in.close();	
		}
		catch(IOException e){
		  e.printStackTrace();
		}
	 }
	 public void printAllSequences(){
         for(SequenceFasta fasta: sequenceList){
             DnaSequenceFasta fastaDna = (DnaSequenceFasta) fasta;
			 System.out.print(fastaDna);
		 }
	 }

	 public void printSequenceStatistics(){
		 int minLen = Integer.MAX_VALUE, maxLen = 0;
		 int avgLen = 0;
		 int len = 0, sum = 0;
		 int size = sequenceList.size();
		 ArrayList<Integer> intList = new ArrayList<Integer>();
		  
		 for(SequenceFasta fasta: sequenceList){
            DnaSequenceFasta fastaDna = (DnaSequenceFasta) fasta;
		    len = fastaDna.getLength();
			intList.add(Integer.valueOf(len));
            if(len < minLen) minLen = len;
			if(len > maxLen) maxLen = len;
			sum += len;
		 }
		 Collections.sort(intList);
         int i = 0;
		 int sum10 = 0, sum100 = 0, N50 = 0, tmp = 0;
		 int cutoff = (int) (sum * .5);
		 boolean flag = false;
		 for(Integer length: intList){
            ++i;
			if(i <= 10) sum10 += length;
			if(i <=100) sum100 += length;
			tmp += length;
			if(!flag && tmp > cutoff){
				N50 = len;
				flag = true;
			}
			if(i > 100 && flag) break;
		 }
		 avgLen = Math.round( (float) sum / size);
		 System.out.printf("There are %d sequences.\n", size);
		 System.out.printf("Total length is %d\n", sum);
		 System.out.printf("Minimum length of sequences is %d.\n", minLen);
		 System.out.printf("Mamimum length of sequences is %d.\n", maxLen);
		 System.out.printf("Average length of sequences is %d.\n", avgLen);
		 System.out.printf("The length of largest 10 sequences is %d.\n", sum10);
		 System.out.printf("The length of largest 100 sequences is %d.\n", sum100);
		 System.out.printf("N50 is %d.\n", N50);
	 }
	 public void sortSequenceByLen(){
         Collections.sort(sequenceList);
	 }
	 public void reverseSequence(){
         Collections.reverse(sequenceList);
	 }
	 //For testing only
	 public static void main(String [] args){
		   SequenceList list = new SequenceList(); 
           list.getFromFileFasta("test2.fna");
		   //list.sortSequenceByLen();
		   //list.reverseSequence();
           //list.printAllSequences();
		   list.printSequenceStatistics();
	 }
}
