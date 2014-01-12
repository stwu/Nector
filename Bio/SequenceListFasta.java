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
 * An object of <code>SequenceListFasta</code> is a 
 * collectin of dna sequences in FASTA format
 * @author Sitao Wu
 * @version 0.1
 * @since  2013
 */
public class SequenceListFasta{
     private LinkedList<SequenceFasta> sequenceList;//Sequence List
	 private String inputFile;//input file name
	 private Scanner scanner;//input stream
	 private String oldLine;//store previous line
	 public final static int CAPACITY = 400;//used for StringBuilder
	 public final static int MAXSEQ = 10000; //cache: read 10,000 sequences at a time
	 
	 /**
	   * SequenceListFasta constructor
	   */
     public SequenceListFasta(){
        sequenceList = new LinkedList<SequenceFasta>();
		inputFile = null;
		scanner = null;
		oldLine = null;
	 }
     
     /**
	   * SequenceListFasta constructor
	   * @param fileName for input file in FASTA format
	   */
	 public SequenceListFasta(String fileName){
		this();
		inputFile = fileName; 
	 }
	 
	 /**
	   * Clear all the members in SequenceListFasta object    
	   */
	 public void clear(){
        sequenceList.clear();       
		scanner = null;
		oldLine = null;
	 }

	 /**
	   * Get the sequence list
	   * @return the list of sequences
	   */
	 public LinkedList<SequenceFasta> getList(){
        return sequenceList;
	 }

	 /**
	   * Get the input file name
	   * @return the file name of input data
	   */
	 public String getFile(){
        return inputFile;
	 }

	 /**
	   * Get the scanner object
	   * @return the scanner object
	   */	 
	 protected Scanner getScanner(){
        return scanner;
	 }

	 /**
	   * Get the last line when reading input file
	   * @return the last line when reading input file
	   */	 
	 protected String getOldLine(){
        return oldLine;
	 }
	 
	 /**
	   * Set the scanner object
	   * @param sc is the scanner object
	   */	  
     protected void setScanner(Scanner sc){
		 scanner = sc;
	 }
     
	 /**
	   * Set the last line when reading input file
	   * @param oLine is the last line when reading input file
	   */	
     protected void setOldLine(String oLine){
		 oldLine = oLine;
	 }
     
	 /**
	   * Excluding/selecting reads in a list file
	   * @param listFile is a list file containing read names to be excluded/selected
	   * @param outputFileName is the output file after excluding/selecting
	   * @param flag is used for: true for excluding; false for selecting 
	   */
     public void excluding(String listFile, String outputFileName, boolean flag){
        HashSet<String> set = new HashSet<String>(); 
		String firstDesc = null;
		String line = null;
		int len;
		clear();
		try{
	      Scanner in = new Scanner(new BufferedReader(new FileReader(listFile)));
	      while(in.hasNextLine()){
		     line = in.nextLine();
			 set.add(line);
          }
		  in.close();
		  boolean finished = false;
		  PrintWriter out = new PrintWriter(outputFileName);
		  while(!finished){
              sequenceList.clear();
              finished = cache();
		      ListIterator<SequenceFasta> it = sequenceList.listIterator();
		      while(it.hasNext()){
                 SequenceFasta fasta = it.next();
			     firstDesc = fasta.getDescription().split("\\s+")[0];
                 len = firstDesc.length();
			     firstDesc = firstDesc.substring(1, len);
			     if(flag){
                    if(set.contains(firstDesc)) it.remove();
			     }
			     else{
                   if(!set.contains(firstDesc)) it.remove();
				 }
			  }//while hasNext
              for(SequenceFasta fasta: sequenceList){
                DnaSequenceFasta fastaDna = (DnaSequenceFasta) fasta;
			    out.print(fastaDna);
		      } 
		  }//while finished
		  out.close();
		}//try
		catch(IOException e){
		  e.printStackTrace();
		}
	 }

	 /**
	   * Reading 10000 sequences at a time. It can be considered as a cache.
	   */     
     public boolean cache(){
		String line, desc = null;
		DnaSequenceFasta fasta = null;
		StringBuilder seq = new StringBuilder(CAPACITY);
		seq.setLength(0);
		int i = 0;
		boolean finished = true;
        try{	   
           if(scanner == null){
	           scanner = new Scanner(new BufferedReader(new FileReader(inputFile)));
	   }
		   if(oldLine != null){
             desc = oldLine;
			 ++i;
		   }
	       while(scanner.hasNextLine()){
		     line = scanner.nextLine();
			 if(line.substring(0, 1).equals(">")){
			 	if(i > 0){
                    fasta = new DnaSequenceFasta(desc, seq.toString());
					sequenceList.add(fasta);
					seq.setLength(0);
					if(i > MAXSEQ){
						finished = false;
						oldLine = line;
					    break;
					}
				}
			    desc = line;//new String(line); 	
				++i;
			 }
			 else{
				seq.append(line);
			 }
		   }//while
		   if(finished){
			 scanner.close();
			 if(i > 0){
                fasta = new DnaSequenceFasta(desc, seq.toString());
			    sequenceList.add(fasta);
		 	    seq.setLength(0);
			 }
		   }
		}
		catch(IOException e){
		  e.printStackTrace();
		}
		return finished;
	 }

	 /**
	   * Get the statistical information for sequences
	   * @param outputFileName is the output file
	   */
	 public void printSequenceStatistics(String outputFileName){
		 int minLen = Integer.MAX_VALUE, maxLen = 0;
		 int avgLen = 0;
		 int len = 0, sum = 0;
		 boolean finished = false;
		 int size = 0;
		 LinkedList<Integer> intList = new LinkedList<Integer>();
		 clear();
         while(!finished){
		   sequenceList.clear();
		   finished = cache();
		   size += sequenceList.size();
		   for(SequenceFasta fasta: sequenceList){
              DnaSequenceFasta fastaDna = (DnaSequenceFasta) fasta;
		      len = fastaDna.getLength();
			  intList.add(Integer.valueOf(len));
              if(len < minLen) minLen = len;
			  if(len > maxLen) maxLen = len;
			  sum += len;
		   }
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

	 /**
	   * Print all sequences in FASTA format for sorting and reversing tools
	   * @param outputFileName is the output file
	   */
     public void printAllSequence(String outputFileName){//used for sortSequenceByLen reverseSequence
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
	 /**
	   * Processing sorting and reversing. It needs reading all seqeunces.
	   * @param cmd is the tool name: sort or reverse
	   * @param args is the argument for the tools
	   */
     private void processAll(String cmd, String [] args){
		 boolean finished = false;
		 clear();
		 sequenceList.clear();
		 String outputFileName = null;
		 String order = null;
		 while(!finished){
            finished = cache();
		 }
		 if(cmd.equals("sort")){
            if(args.length < 2){
                System.err.println("Error! For sorting tool, we need two parameters: outputFileName and order!");
				System.exit(1);
			}
            outputFileName = args[0];
			order = args[1];
		    if(order.equals("a")){
               Collections.sort(sequenceList);
		    }
		    else if(order.equals("d")){
               Collections.sort(sequenceList, Collections.reverseOrder());
		    }
		 }
		 else if(cmd.equals("reverse")){
            if(args.length < 1){
                System.err.println("Error! For reversing tool, we need one parameter: outputFileName!");
				System.exit(1);
			}
            outputFileName = args[0];
            Collections.reverse(sequenceList);
		 }
	     printAllSequence(outputFileName);
	 }

	 /**
	   * Sorting sequences by sequence length
	   * @param outputFileName is the output file
	   * @param order is the sorting order: a:ascending;d:descending
	   */
	 public void sortSequenceByLen(String outputFileName, String order){
         processAll("sort", new String [] {outputFileName, order});
	 }

	 /**
	   * Get the reverse order of a sequence list
	   * @param outputFileName is the output file
	   */
	 public void reverseSequence(String outputFileName){
		 processAll("reverse", new String [] {outputFileName});
	 }
	 
	 /**
	   * Processing using cache. Don't need reading all sequences before processing.
	   * @param cmd is the tool name other than sort or reverse
	   * @param args is the argument for the tools
	   */
     private void processCache(String cmd, String [] args){
		 boolean finished = false;
		 String outputFileName = null;
		 String barCode = null;
		 int minLen = 0;
		 clear();
		 if(cmd.equals("filter")){
             if(args.length < 2){
                 System.err.println("Error! For filtering tool, we need two parameters: outputFileName and minLen!");
		         System.exit(1);
			 }
             outputFileName = args[0];
			 minLen = Integer.parseInt(args[1]);
		 }
		 else if(cmd.equals("barcode")){
             if(args.length < 2){
                 System.err.println("Error! For barcode tool, we need two parameters: outputFileName and barCode!");
		         System.exit(1);
			 }
             outputFileName = args[0];
			 barCode = args[1];
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
              sequenceList.clear();
              finished = cache();
              ListIterator<SequenceFasta> it = sequenceList.listIterator();
		      while(it.hasNext()){
                SequenceFasta fasta = it.next();
		        if(cmd.equals("filter")){
                     if(fasta.getLength() < minLen) it.remove();
			    }
				else if(cmd.equals("cr")){
			         it.set((SequenceFasta) ((DnaSequenceFasta) fasta).complementaryReverse());
				}
				else if(cmd.equals("barcode")){
					 DnaSequenceFasta tmp = ((DnaSequenceFasta) fasta).barcode(barCode);
					 if(tmp == null){
					   it.remove();	 
					 }
					 else{
                       it.set((SequenceFasta) tmp);
					 }
				}
		      }
              for(SequenceFasta fasta: sequenceList){
                 DnaSequenceFasta fastaDna = (DnaSequenceFasta) fasta;
			     out.print(fastaDna);
			  }
		   }//while
		   out.close();
		 }
		 catch(IOException e){
			 e.printStackTrace();
		 }

	 }
     
	 /**
	   * Get the complementary reverse of a sequence list
	   * @param outputFileName is the output file
	   */     
	 public void complementaryReverseSequence(String outputFileName){
		 processCache("cr", new String [] {outputFileName});
	 }

	 /**
	   * Filtering short sequence by sequence length
	   * @param outputFileName is the output file
	   * @param minLen is the minimum sequence length
	   */
	 public void filterSequence(String outputFileName, int minLen){
		 processCache("filter", new String [] {outputFileName, String.valueOf(minLen)});
	 }
	 
	 /**
	   * print sequences using cache
	   * @param outputFileName is the output file
	   */	 
	 public void printAllSequenceWithCache(String outputFileName){
		 processCache("print", new String [] {outputFileName});
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
	   * @return true if FASTA format; false otherwise
	   */		 
	 public boolean formatCheck(){
		 boolean finished = false;
		 clear();
		 while(!finished){
            sequenceList.clear();
            finished = cache();
            ListIterator<SequenceFasta> it = sequenceList.listIterator();
		    while(it.hasNext()){
              SequenceFasta fasta = it.next();
		      if(!(((DnaSequenceFasta)fasta).formatCheck())){
                 return false;
			  }
		    }   
		 }
		 return true;
	 }
	 public static void main(String [] args){
		SequenceListFasta list = new SequenceListFasta("example.fasta");
        list.printSequenceStatistics("example.stat.txt");
		list.reverseSequence("example.revere.txt");
		list.complementaryReverseSequence("example.cr.txt");
		list.sortSequenceByLen("ascending", "example.sort.txt");
		list.excluding("list", "example.exclude.txt", true);
		list.filterSequence("example.filter.txt", 100);
		System.out.println(list.formatCheck());
	 }
}
