/* The MIT License

Copyright (c) 2013, by Sitao Wu <wusitao2000@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

package Bio;
import java.io.*;
import java.util.*;
import java.util.LinkedList;
import java.util.Collections;
import java.util.List;

/**
 * An object of <code>SequenceListFastq</code> is a 
 * collectin of dna sequences in FASTQ format
 * @author Sitao Wu
 * @version 0.1
 * @since  2013
 */
public class SequenceListFastq extends SequenceListFasta{

	 /**
	   * SequenceListFastq constructor
	   */
     public SequenceListFastq(){
		super();
	 }

     /**
	   * SequenceListFastq constructor
	   * @param fileName for input file in FASTQ format
	   */
	 public SequenceListFastq(String fileName){
		super(fileName);
	 }
	 
	 /**
	   * Reading 10000 sequences at a time. It can be considered as a cache.
	   */   
	 @Override
	 public boolean cache(){
		String line, desc = null,qDesc = null;
		DnaSequenceFastq fastq = null;
		StringBuilder seq = new StringBuilder(CAPACITY);
		StringBuilder qSeq = new StringBuilder(CAPACITY);
		seq.setLength(0);
		qSeq.setLength(0);
		int i = 0;
		boolean finished = true;
		byte flagSeq = 0;
		String oldLine = super.getOldLine();
		Scanner scanner = getScanner();
		LinkedList<SequenceFasta> sequenceList = super.getList();
        try{
           if(scanner == null){
	           super.setScanner(new Scanner(new BufferedReader(new FileReader(super.getFile()))));
			   scanner = super.getScanner();
		   }
		   if(oldLine != null){
			 flagSeq = 1;
             desc = oldLine;
			 //System.out.println("desc=" + desc);
			 ++i;
		   }
	       while(scanner.hasNextLine()){
		     line = scanner.nextLine();
		     //System.out.println(line);
			 if((flagSeq == 0) && line.substring(0, 1).equals("@")){
			 	if(i > 0){
                    fastq = new DnaSequenceFastq(desc, seq.toString(), qDesc, qSeq.toString());
					sequenceList.add(fastq);
					seq.setLength(0);
					qSeq.setLength(0);
					if(i > MAXSEQ){
						finished = false;
			            //System.out.println("desc=" + line);
						super.setOldLine(line);
					    break;
					}
				}
				flagSeq = 1;
			    desc = line; 	
				++i;
			 }
			 else if( (flagSeq == 1) && line.substring(0, 1).equals("+")){
				flagSeq = 2;
                qDesc = line;
			 }
			 else{
                if(flagSeq == 1){
				   seq.append(line);
				}
				else if(flagSeq == 2){
				   qSeq.append(line);
				   if(qSeq.length() == seq.length()){
					   flagSeq = 0;					   
				   }
				}
			 }
		   }//while
		   if(finished){
			 scanner.close();
			 if(i > 0){
                fastq = new DnaSequenceFastq(desc, seq.toString(), qDesc, qSeq.toString());
			    sequenceList.add(fastq);
		 	    seq.setLength(0);
		 	    qSeq.setLength(0);
			 }
		   }
		}
		catch(IOException e){
		  e.printStackTrace();
		}
		return finished;
	 }

	 /**
	   * Print all sequences in FASTQ format for sorting and reversing tools
	   * @param outputFileName is the output file
	   */
	 @Override
     public void printAllSequence(String outputFileName){//used for sortSequenceByLen reverseSequence
	     LinkedList<SequenceFasta> list = null;
		 try{
		   PrintWriter out = new PrintWriter(outputFileName);
		   list = super.getList();
           for(SequenceFasta fasta: list){
              DnaSequenceFastq fastqDna = (DnaSequenceFastq) fasta;
			  out.print(fastqDna);
		   }
		   out.close();
		 }
		 catch(IOException e){
			 e.printStackTrace();
		 }
	 }

	 /**
	   * Processing using cache. Don't need reading all sequences before processing.
	   * @param cmd is the tool name other than sort or reverse
	   * @param args is the argument for the tools
	   */	 
     private void processCache(String cmd, String [] args){
		 boolean finished = false;
		 LinkedList<SequenceFasta> list = super.getList();
		 super.clear();
		 String outputFileName = null;
		 String barCode = null;
		 int cutoff = 0;
		 int minLen = 0;
		 float maxPer = 0.0F;
		 if(cmd.equals("trim")){
             if(args.length < 4){
                 System.err.println("Error! For trimming tool, we need four parameters: outputFileName, cutoff, maxPer and minLen!");
		         System.exit(1);
			 }
			 outputFileName = args[0];
			 cutoff = Integer.parseInt(args[1]);
			 maxPer = Float.parseFloat(args[2]);
			 minLen = Integer.parseInt(args[3]);
		 }
		 else if(cmd.equals("barcode")){
             if(args.length < 2){
                 System.err.println("Error! For barcode tool, we need two parameters: outputFileName and barcode!");
		         System.exit(1);
			 }
			 outputFileName = args[0];
			 barCode = args[1];
		 }
		 else if(cmd.equals("convert")){
             if(args.length < 1){
                 System.err.println("Error! For converting tool, we need one parameter: outputFileName!");
		         System.exit(1);
			 }
			 outputFileName = args[0];
		 }
		 else if(cmd.equals("cr")){
             if(args.length < 1){
                 System.err.println("Error! For complementaryReverse tool, we need one parameter: outputFileName!");
		         System.exit(1);
			 }
			 outputFileName = args[0];
		 }
		 else if(cmd.equals("print")){
             if(args.length < 1){
                 System.err.println("Error! For print tool, we need one parameter: outputFileName!");
		         System.exit(1);
			 }
			 outputFileName = args[0];
		 }
		 try{
		   PrintWriter out = new PrintWriter(outputFileName);
		   while(!finished){
              list.clear();
			  //System.out.println("-------");
              finished = cache();
			  //System.out.println("=======" + finished);
              ListIterator<SequenceFasta> it = list.listIterator();
		      while(it.hasNext()){
                SequenceFasta fasta = it.next();
				if(cmd.equals("trim")){
                   ((DnaSequenceFastq)fasta).trim(cutoff, maxPer, minLen);
			       if(fasta.getLength() == 0) it.remove();
				}
				else if(cmd.equals("barcode")){
				   DnaSequenceFastq tmp = ((DnaSequenceFastq) fasta).barcode(barCode);
				   if(tmp == null){
					  it.remove();	 
				   }
				   else{
                      it.set((SequenceFasta) tmp);
				   }
				}
				else if(cmd.equals("convert")){
			       it.set((SequenceFasta) ((DnaSequenceFastq) fasta).convertToFasta());
				}
				else if(cmd.equals("cr")){
			       it.set((SequenceFasta) ((DnaSequenceFastq) fasta).complementaryReverse());
				}
		      }   
              for(SequenceFasta fasta: list){
				 if(cmd.equals("convert")){
                   DnaSequenceFasta fastaDna = (DnaSequenceFasta) fasta;
			       out.print(fastaDna);
				 }
				 else{
                   DnaSequenceFastq fastqDna = (DnaSequenceFastq) fasta;
			       out.print(fastqDna);
				 }
			  }
			  //out.println("-------");
			  //System.out.println("0000000000");
		   }
		   out.close();
		 }
		 catch(IOException e){
			 e.printStackTrace();
		 }
      
	 }
	 /**
	  * Triming reads with low quality tails for illumina reads
	  * @param cutoff is the minimum quality score cutoff(default value = 13)
	  * @param maxPer is the maximum percentage of mismatches
	  * @param minLen is the minimum length after trimming
	  */
	 public void trim(String outputFileName, int cutoff, float maxPer, int minLen){
         processCache("trim", new String [] {outputFileName, String.valueOf(cutoff), String.valueOf(maxPer), String.valueOf(minLen)});
	 }
	 /**
	  * Converting FASTQ format to FASTA format
	  */
	 public void convertToFasta(String outputFileName){
         processCache("convert", new String [] {outputFileName});
	 }
	 /**
	   * print sequences using cache
	   * @param outputFileName is the output file
	   */	
	 public void printAllSequenceWithCache(String outputFileName){
		 processCache("print", new String [] {outputFileName});
	 }
	 /**
	   * Get the complementary reverse of a sequence list
	   * @param outputFileName is the output file
	   */	 
	 @Override
	 public void complementaryReverseSequence(String outputFileName){
         processCache("cr", new String [] {outputFileName});
	 }
	 /**
	   * Extracting sequence with bar code using cache
	   * @param outputFileName is the output file
	   */		 	 
	 public void barCode(String outputFileName, String barcode){
		 processCache("barcode", new String [] {outputFileName, barcode});
	 }
	 /**
	   * Checking the format 
	   * @return true if FASTQ format; false otherwise
	   */		
	 @Override
	 public boolean formatCheck(){
		 boolean finished = false;
		 LinkedList<SequenceFasta> list = super.getList();
		 super.clear();
		 while(!finished){
            list.clear();
            finished = cache();
            ListIterator<SequenceFasta> it = list.listIterator();
		    while(it.hasNext()){
              SequenceFasta fasta = it.next();
		      if(!(((DnaSequenceFastq)fasta).formatCheck())){
                 return false;
			  }
		    }   
		 }
		 return true;
	 }
	 public static void main(String [] args){
		SequenceListFastq list = new SequenceListFastq("example.fastq");
        list.printSequenceStatistics("example.stat.txt");
		list.reverseSequence("example.revere.txt");
		list.complementaryReverseSequence("example.cr.txt");
		list.sortSequenceByLen("ascending", "example.sort.txt");
		list.excluding("list", "example.exclude.txt", true);
		list.filterSequence("example.filter.txt", 100);
		list.convertToFasta("example.fasta.txt");
		list.trim("example.trim.txt", 13, .05F, 50);
		System.out.println(list.formatCheck());
	 }
}
