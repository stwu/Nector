package Bio;
import java.io.*;
import java.util.*;
import java.util.LinkedList;
import java.util.Collections;
import java.util.List;

public class SequenceListFastq extends SequenceListFasta{

     public SequenceListFastq(){
		super();
	 }
	 public SequenceListFastq(String fileName){
		this();
		getFromFile(fileName); 
	 }
	 public void getFromFile(String fileName){
		String line, line2, line3, line4, desc = null;
		DnaSequenceFastq fastq = null;
		int i = 0;
		try{
	      Scanner in = new Scanner(new BufferedReader(new FileReader(fileName)));
	      while(in.hasNextLine()){
		     line = in.nextLine();
		     line2 = in.nextLine();
		     line3 = in.nextLine();
		     line4 = in.nextLine();
			 if(line.substring(0, 1).equals("@")){
				super.getList().add(new DnaSequenceFastq(line, line2, line3, line4));
			 }
		  }
	      in.close();	
		}
		catch(IOException e){
		  e.printStackTrace();
		}
	 }



	 public void printAllSequencesFasta(String outputFileName){
		 LinkedList<SequenceFasta> list = super.getList(); 
		 try{
		   PrintWriter out = new PrintWriter(outputFileName);
		   list = super.getList();
           for(SequenceFasta fasta: list){
			   out.print(fasta);
		   }
		   out.close();
		 }
		 catch(IOException e){
			 e.printStackTrace();
		 }
     }
	 public void printAllSequences(String outputFileName){
		 DnaSequenceFastq fastqDna = null;
		 LinkedList<SequenceFasta> list = super.getList(); 
		 try{
		   PrintWriter out = new PrintWriter(outputFileName);
		   list = super.getList();
           for(SequenceFasta fasta: list){
               fastqDna = (DnaSequenceFastq) fasta;
			   out.print(fastqDna);
		   }
		   out.close();
		 }
		 catch(IOException e){
			 e.printStackTrace();
		 }
	 }

	 public void trim(int cutoff, float maxPer, int minLen){
         LinkedList<SequenceFasta> list = super.getList();
		 ListIterator<SequenceFasta> it = list.listIterator();
		 SequenceFasta fasta = null;
		 while(it.hasNext()){
             fasta = it.next();
             ((DnaSequenceFastq)fasta).trim(cutoff, maxPer, minLen);
			 if(fasta.getLength() == 0) it.remove();
		 }
	 }
	 public boolean formatCheck(){
		 ListIterator<SequenceFasta> it = super.getList().listIterator();
		 SequenceFasta fasta = null;
		 while(it.hasNext()){
              fasta = it.next();
			  if(!(((DnaSequenceFastq)fasta).formatCheck())){
                 return false;
			  }
		 }
		 return true;
	 }
	 public void convertToFasta(){
         LinkedList<SequenceFasta> list = super.getList();
		 ListIterator<SequenceFasta> it = list.listIterator();
		 SequenceFasta fasta = null;
		 while(it.hasNext()){
             fasta = it.next();
			 it.set((SequenceFasta) ((DnaSequenceFastq) fasta).convertToFasta()); 
		 }
	 }
	 public void complementaryReverseSequence(){
		 LinkedList<SequenceFasta> list = super.getList();
		 ListIterator<SequenceFasta> it = list.listIterator();
		 SequenceFasta fasta = null;
		 while(it.hasNext()){
              fasta = it.next();
			  it.set((SequenceFasta) ((DnaSequenceFastq) fasta).complementaryReverse());
		 }
	 }
}
