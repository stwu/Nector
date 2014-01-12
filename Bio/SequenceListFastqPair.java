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

/**
 * An object of <code>SequenceListFastqPair</code> is a 
 * collectin of dna pair-ended sequences in FASTQ format
 * @author Sitao Wu
 * @version 0.1
 * @since  2013
 */

public class SequenceListFastqPair extends SequenceListFastaPair{
    //private SequenceListFasta [] pair = new SequenceListFasta[2];
	
	 /**
	   * SequenceListFastqPair constructor
	   */
    public SequenceListFastqPair(){
		super();
       //pair[0] = null;
       //pair[1] = null;
	}
	 /**
	   * SequenceListFastaPair constructor
	   * @param pairFileName for pair-ended input files in FASTQ format
	   */   
	public SequenceListFastqPair(String pairFileName){
	   //System.out.println(pairFileName);
       String [] files = pairFileName.split(",");
	   if(files.length != 2){
		   System.err.println("Error! Sequence pair shoulde have two input files!");
		   System.exit(1);   
	   }
	   super.setPair(new SequenceListFastq(files[0]), new SequenceListFastq(files[1]));
	   //System.out.println(files[0]);
	   //System.out.println(files[1]);
       //pair[0] = new SequenceListFastq(files[0]);
       //pair[1] = new SequenceListFastq(files[1]);
	}
	/**
	 * Assembly pair-ended reads into  longer reads
	 * @param maxMismatch is the maximum number of mismatches
	 * @param minOverlap is the minimum lemgth of overlapped region
	 */
	@Override
	public void assembly(String outputFileName, int maxMismatch, int minOverlap){
		boolean finished1 = false;
		boolean finished2 = false;
		try{
		  PrintWriter out = new PrintWriter(outputFileName);
		  while(!finished1 && !finished2){
              SequenceListFasta [] myPair = super.getPair();
              LinkedList<SequenceFasta> list1 = myPair[0].getList();
              LinkedList<SequenceFasta> list2 = myPair[1].getList();
			  list1.clear();
			  list2.clear();
              finished1 = myPair[0].cache();
              finished2 = myPair[1].cache();
		      ListIterator<SequenceFasta> it1 = list1.listIterator();
		      ListIterator<SequenceFasta> it2 = list2.listIterator();
		      while(it1.hasNext() && it2.hasNext()){
				 DnaSequenceFastq fastq1 = (DnaSequenceFastq) (it1.next()); 
				 DnaSequenceFastq fastq2 = (DnaSequenceFastq) (it2.next()); 
                 DnaSequenceFasta fasta1 = fastq1.convertToFasta();
                 DnaSequenceFasta fasta2 = fastq2.convertToFasta();
				 //System.out.println(fastq1); 
				 if(!fasta1.assembly(fasta2, maxMismatch, minOverlap)){//does not assembly
                   it1.remove();
				 }
				 else{
                   it1.set(fasta1);
				 }
			  }//while hasNext
              for(SequenceFasta fasta1: list1){
                DnaSequenceFasta fastaDna = (DnaSequenceFasta) fasta1;
			    out.print(fastaDna);
		      }
          }//while
		  out.close();
	  }
	  catch(IOException e){
		 e.printStackTrace();
      }
   }
   public static void main(String [] args){
      SequenceListFastqPair p = new SequenceListFastqPair("read1.fq,read2.fq");
	  p.assembly("ooo2", 3, 10);
   }
}
