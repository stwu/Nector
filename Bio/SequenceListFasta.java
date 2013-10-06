package Bio;
import java.io.*;
import java.util.*;
import java.util.LinkedList;
import java.util.Collections;
import java.util.List;

public class SequenceListFasta{
     private LinkedList<SequenceFasta> sequenceList;
	 public final static int CAPACITY = 400;

     public SequenceListFasta(){
        sequenceList = new LinkedList<SequenceFasta>();
	 }
	 public SequenceListFasta(String fileName){
		this();
		getFromFile(fileName); 
	 }
	 public LinkedList<SequenceFasta> getList(){
        return sequenceList;
	 }
     public void excluding(String listFile, boolean flag){
        HashSet<String> set = new HashSet<String>(); 
		String firstDesc = null;
		String line = null;
		int len;
		try{
	      Scanner in = new Scanner(new BufferedReader(new FileReader(listFile)));
	      while(in.hasNextLine()){
		     line = in.nextLine();
			 set.add(line);
          }
		  in.close();
		  ListIterator<SequenceFasta> it = sequenceList.listIterator();
		  SequenceFasta fasta = null;
		  while(it.hasNext()){
              fasta = it.next();
			  firstDesc = fasta.getDescription().split("\\s+")[0];
              len = firstDesc.length();
			  firstDesc = firstDesc.substring(1, len);
			  if(flag){
                 if(set.contains(firstDesc)) it.remove();
			  }
			  else{
                 if(!set.contains(firstDesc)) it.remove();
			  }
		  }
		}
		catch(IOException e){
		  e.printStackTrace();
		}
		catch(Exception e){
		  e.printStackTrace();
		}

	 }
	 public void getFromFile(String fileName){
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
	 public void printAllSequences(String outputFileName){
		 try{
		   PrintWriter out = new PrintWriter(outputFileName);
           for(SequenceFasta fasta: sequenceList){
               DnaSequenceFasta fastaDna = (DnaSequenceFasta) fasta;
			   out.print(fastaDna);

		   }
		   out.close();
		 }
		 catch(IOException e){
			 e.printStackTrace();
		 }
	 }

	 public void printSequenceStatistics(String outputFileName){
		 int minLen = Integer.MAX_VALUE, maxLen = 0;
		 int avgLen = 0;
		 int len = 0, sum = 0;
		 int size = sequenceList.size();
		 LinkedList<Integer> intList = new LinkedList<Integer>();
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
		 try{
		    PrintWriter out = new PrintWriter(outputFileName);
		    out.printf("There are %d sequences.\n", size);
		    out.printf("Total length is %d\n", sum);
		    out.printf("Minimum length of sequences is %d.\n", minLen);
		    out.printf("Mamimum length of sequences is %d.\n", maxLen);
		    out.printf("Average length of sequences is %d.\n", avgLen);
		    out.printf("The length of largest 10 sequences is %d.\n", sum10);
		    out.printf("The length of largest 100 sequences is %d.\n", sum100);
		    out.printf("N50 is %d.\n", N50);
			out.close();
		 }
		 catch(IOException e){
            e.printStackTrace();
		 }
	 }
     
	 public void sortSequenceByLen(String order){
		 if(order.equals("ascending")){
            Collections.sort(sequenceList);
		 }
		 else if(order.equals("descending")){
            Collections.sort(sequenceList, Collections.reverseOrder());
		 }
	 }
	 public void reverseSequence(){
         Collections.reverse(sequenceList);
	 }
	 public void complementaryReverseSequence(){
		 ListIterator<SequenceFasta> it = sequenceList.listIterator();
		 SequenceFasta fasta = null;
		 while(it.hasNext()){
              fasta = it.next();
			  it.set((SequenceFasta) ((DnaSequenceFasta) fasta).complementaryReverse());
		 }
	 }
	 public void filterSequence(int minLen){
		 ListIterator<SequenceFasta> it = sequenceList.listIterator();
		 SequenceFasta fasta = null;
		 while(it.hasNext()){
              fasta = it.next();
              if(fasta.getLength() < minLen) it.remove();
		 }
	 }
	 public boolean formatCheck(){
		 ListIterator<SequenceFasta> it = sequenceList.listIterator();
		 SequenceFasta fasta = null;
		 while(it.hasNext()){
              fasta = it.next();
			  if(!(((DnaSequenceFasta)fasta).formatCheck())){
                 return false;
			  }
		 }
		 return true;
	 }
}
