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
 * An object of <code>SequenceListFastaPair</code> is a 
 * collectin of dna pair-ended sequences in FASTA format
 * @author Sitao Wu
 * @version 0.1
 * @since  2013
 */
public class SequenceListFastaPair{
    private SequenceListFasta [] pair = new SequenceListFasta[2];
	 /**
	  * SequenceListFastaPair constructor
	  */
    public SequenceListFastaPair(){
       pair[0] = null;
       pair[1] = null;
	}
	 /**
	  * Set the pair-ened file names
	  * @param f1 for forward input file in FASTA format
	  * @param f2 for reverse input file in FASTA format
	  */    
	public void setPair(SequenceListFasta f1, SequenceListFasta f2){
       pair[0] = f1;
       pair[1] = f2;
	}
	/**
	 * Get the pair-ended file names
	 * @return the pair-ended file names in an array
	 */
	public SequenceListFasta [] getPair(){
       return pair;
	}
	/**
	 * SequenceListFastaPair constructor
	 * @param pairFileName for pair-ended input files in FASTA format
	 */
	public SequenceListFastaPair(String pairFileName){
       String [] files = pairFileName.split(",");
	   if(files.length != 2){
		   System.err.println("Error! Sequence pair shoulde have two input files!");
		   System.exit(1);   
	   }
       pair[0] = new SequenceListFasta(files[0]);
       pair[1] = new SequenceListFasta(files[1]);
	}
	/**
	 * Assembly pair-ended reads into  longer reads
	 * @param maxMismatch is the maximum number of mismatches
	 * @param minOverlap is the minimum lemgth of overlapped region
	 */
	public void assembly(String outputFileName, int maxMismatch, int minOverlap){
		boolean finished1 = false;
		boolean finished2 = false;
		try{
		  PrintWriter out = new PrintWriter(outputFileName);
		  while(!finished1 && !finished2){
              LinkedList<SequenceFasta> list1 = pair[0].getList();
              LinkedList<SequenceFasta> list2 = pair[1].getList();
			  list1.clear();
			  list2.clear();
              finished1 = pair[0].cache();
              finished2 = pair[1].cache();
		      ListIterator<SequenceFasta> it1 = list1.listIterator();
		      ListIterator<SequenceFasta> it2 = list2.listIterator();
		      while(it1.hasNext() && it2.hasNext()){
                 SequenceFasta fasta1 = it1.next();
                 SequenceFasta fasta2 = it2.next();
				 
				 if(!((DnaSequenceFasta) fasta1).assembly( (DnaSequenceFasta) fasta2, maxMismatch, minOverlap)){//does not assembly
                   it1.remove();
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
      SequenceListFastaPair p = new SequenceListFastaPair("read1.fa,read2.fa");
	  p.assembly("ooo", 3, 10);
   }
}
