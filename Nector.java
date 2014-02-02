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

import Bio.*;
import java.io.*;
import java.util.*;

/**
 * <code>Nector</code> is a class for  
 * sequence preprocessing tool on NGS data
 * @author Sitao Wu
 * @version 0.1
 * @since 2013
 */
public class Nector{
    private enum Format{
		 OTHER,
		 FASTA,
		 FASTQ
    };
    private static Format format;
	private static String version = "0.1";
	/**
	 * Main program starts here
         * @param args for program arguments 
	 */
    public static void main(String [] args){
		processArg(args);
    }
	/**
	 * Providing help information for running the main program
	 */
	public static void help(){
		System.out.println("");
		System.out.println("Program: Nector - sequence preprocessing tools for NGS data");
		System.out.printf("Version %s. Developed by Sitao Wu\n", version);
		System.out.println("The MIT license");
		System.out.println("");
        System.out.println("Usage:\tjava -jar Nector.jar toolName inputFileName(s) outputFileName parameters");
		System.out.println("");
		System.out.println("Tools supported:");
		System.out.println("");
		System.out.println("1.stat: Get statistics of sequence file");
		System.out.println("");
		System.out.println("\tExample 1: java -jar Nector.jar stat example.fasta stat.output");
		System.out.println("\tExample 2: java -jar Nector.jar stat example.fastq stat.output");
		System.out.println("");
		System.out.println("2.rev: Get sequences in revserse order");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("");
		System.out.println("\tExample 1: java -jar Nector.jar rev example.fasta exampleReverse.fasta [-l 50]");
		System.out.println("\tExample 2: java -jar Nector.jar rev example.fastq exampleReverse.fastq [-l 50]");
		System.out.println("");
		System.out.println("3.cr: Get sequences in complementary reverse");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("");
		System.out.println("\tExample 1: java -jar Nector.jar cr example.fasta exampleCr.fasta [-l 50]");
		System.out.println("\tExample 2: java -jar Nector.jar cr example.fastq exampleCr.fastq [-l 50]");
		System.out.println("");
		System.out.println("4.sort: Get sequences in sorted order (ascending or descending) by length");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("\tParameter: -o a:ascending;d:descending (default=d)");
		System.out.println("");
		System.out.println("\tExample 1: java -jar Nector.jar sort example.fasta exampleSort.fasta [-l 50] (default sort order = descending)");
		System.out.println("\tExample 2: java -jar Nector.jar sort example.fastq exampleSort.fastq [-l 50] (default sort order = descending)");
		System.out.println("\tExample 3: java -jar Nector.jar sort example.fasta exampleSort.fasta [-o a] [-l 50]");
		System.out.println("\tExample 4: java -jar Nector.jar sort example.fastq exampleSort.fastq [-o a] [-l 50]");
		System.out.println("\tExample 5: java -jar Nector.jar sort example.fasta exampleSort.fasta [-o d] [-l 50]");
		System.out.println("\tExample 6: java -jar Nector.jar sort example.fastq exampleSort.fastq [-o d] [-l 50]");
		System.out.println("");
		System.out.println("5.filter: filter short sequences");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("\tParameter: -n minimum length (default=40)");
		System.out.println("");
		System.out.println("\tExample 1: java -jar Nector.jar filter example.fasta exampleFiltered.fasta [-n 100] [-l 50] (minimum sequence length=100)");
		System.out.println("\tExample 2: java -jar Nector.jar filter example.fastq exampleFiltered.fasta [-n 100] [-l 50] (minimum sequence length=100)");
		System.out.println("");
		System.out.println("6.exclude: excluding sequences in a list");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("");
		System.out.println("\tExample 1: java -jar Nector.jar exclude example.fasta,list exampleExcluded.fasta [-l 50] (excluding reads defined in list by read name)");
		System.out.println("\tExample 2: java -jar Nector.jar exclude example.fastq,list exampleExcluded.fastq [-l 50] (excluding reads defined in list by read name)");
		System.out.println("");
		System.out.println("7.select: selecting sequences in a list");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("");
		System.out.println("\tExample 1: java -jar Nector.jar select example.fasta,list exampleSelecteded.fasta [-l 50] (selecting reads defined in list by read name)");
		System.out.println("\tExample 2: java -jar Nector.jar select example.fastq,list exampleSelecteded.fastq [-l 50] (selecting reads defined in list by read name)");
		System.out.println("");
		System.out.println("8.fastq2fasta: converting fastq to fasta format");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("");
		System.out.println("\tExample: java -jar Nector.jar fastq2fasta example.fastq exampleConvert.fasta [-l 50]");
		System.out.println("");
		System.out.println("9.check: check format");
		System.out.println("");
		System.out.println("\tExample 1: java -jar Nector.jar check example.fasta");
		System.out.println("\tExample 2: java -jar Nector.jar check example.fastq");
		System.out.println("");
		System.out.println("10.trim: trimming and filtering");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("\tParameter: -c minimum quality score (default=13)");
		System.out.println("\tParameter: -p maximum error rate (default=0.05)");
		System.out.println("\tParameter: -L minimum sequence length (default=40");
		System.out.println("");
		System.out.println("\tExample: java -jar Nector.jar trim example.fastq exampleTrimed.fastq [-c 13] [-p 0.05] [-L 70] [-l 50]");
		System.out.println("");
		System.out.println("11.print: write file according to specified lengh per line");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("");
		System.out.println("\tExample: java -jar Nector.jar print example.fasta exampleNew.fasta [-l 50]");
		System.out.println("\tExample: java -jar Nector.jar print example.fastq exampleNew.fastq [-l 50]");
		System.out.println("");
		System.out.println("12.barcode: extract sequences according to barcode");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("\tParameter: -b bar code to be filtered");
		System.out.println("");
		System.out.println("\tExample: java -jar Nector.jar barcode example.fasta exampleBarcode.fasta -b ACC [-l 50]");
		System.out.println("\tExample: java -jar Nector.jar barcode example.fastq exampleBarcode.fastq -b ACC [-l 50]");
		System.out.println("");
		System.out.println("13.assembly: assembly pair-ended reads into linked ones");
		System.out.println("");
		System.out.println("\tParameter: -l line width for sequence (default=250)");
		System.out.println("\tParameter: -m maximum number of mismatch (default=3)");
		System.out.println("\tParameter: -v minimum overlapping length (default=10)");
		System.out.println("");
		System.out.println("\tExample: java -jar Nector.jar assembly example1.fasta,example2.fasta exampleAssembled.fasta [-m 3] [-v 10] [-l 50]");
		System.out.println("\tExample: java -jar Nector.jar assembly example1.fastq,example2.fastq exampleAssembled.fasta [-m 3] [-v 10] [-l 50]");
		System.out.println("");
	}
	
	/**
	 * Run different tools according to arguments 
	 * @param args for program arguments
	 */
	public static void processArg(String [] args){
		String inputFileName = null, listFile = null;
		String outputFileName = null;
		String inputFileName2 = null;
		int n = 0;
        if(args.length < 2){//toolName + inputFileName
		  help();
		  System.exit(1);
		}
        else if(args.length == 2 && !(args[0].equals("check"))){//check tool does not need outputFileName
		  help();
		  System.exit(1);
		}
		inputFileName = args[1];
		if(!(args[0].equals("check"))){
			outputFileName = args[2];
		}
		if(args[0].equals("exclude") || args[0].equals("select")){
		   String [] tmpStr = inputFileName.split(",");
           inputFileName = new String(tmpStr[0]);
           listFile = new String(tmpStr[1]);
		}
		else if(args[0].equals("assembly")){
		   String [] tmpStr = inputFileName.split(",");
           inputFileName = new String(tmpStr[0]);
           inputFileName2 = new String(tmpStr[1]);
		}
		format = fastCheck(inputFileName);
		if(format == Format.OTHER){
			 System.err.println("Error! Input file is not FASTA or FASTQ file!");
             System.exit(1);
		}

        int formatLen = 250; //Default value of line width for sequence

		if(args[0].equals("stat")){//statistics
			runStatistics(inputFileName, outputFileName);
			System.out.println("Nector finished the \"stat\" job successfully!");
		}
		else if(args[0].equals("assembly")){//assembly pair-ended reads into linked ones
			n = 3;
			int maxMismatch = 3;//default value
			int minOverlap = 10;//default value
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else if(args[i].equals("-m")){
                     maxMismatch = Integer.parseInt(args[i + 1]);
				  }
				  else if(args[i].equals("-v")){
                     minOverlap = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
              runAssembly(inputFileName, inputFileName2, outputFileName, maxMismatch, minOverlap);
              System.out.println("Nector finished the \"assembly\" job successfully!");
			}
		}
		else if(args[0].equals("print")){//print FASTA with a specified line length for sequence
			n = 3;
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
              runPrint(inputFileName, outputFileName);
              System.out.println("Nector finished the \"print\" job successfully!");
			}
		}
		else if(args[0].equals("check")){//chekc for FASTA or FASTQ format
            runCheck(inputFileName);
            System.out.println("Nector finished the \"check\" job successfully!");
		}
		else if(args[0].equals("rev")){
			n = 3;
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
			  runReverseOrder(inputFileName, outputFileName);
			  System.out.println("Nector finished the \"rev\" job successfully!");
			}
		}
		else if(args[0].equals("cr")){
			n = 3;
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
			  runComplementaryReverse(inputFileName, outputFileName);
			  System.out.println("Nector finished the \"cr\" job successfully!");

			}
		}
		else if(args[0].equals("trim")){
			n = 3;
			int cutoff = 13;//default value
			float maxPer = 0.05F;//default value
			int minLen = 40;//default value
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else if(args[i].equals("-c")){
                     cutoff = Integer.parseInt(args[i + 1]);
				  }
				  else if(args[i].equals("-p")){
                     maxPer = Float.parseFloat(args[i + 1]);
				  }
				  else if(args[i].equals("-L")){
                     minLen = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
			  runTrim(inputFileName, outputFileName, cutoff, maxPer, minLen);
			  System.out.println("Nector finished the \"trim\" job successfully!");
			}
		}
		else if(args[0].equals("filter")){
			n = 3;
			int minLen = 40;//default value
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else if(args[i].equals("-n")){
                     minLen = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
			  runFilter(inputFileName, outputFileName, minLen);
			  System.out.println("Nector finished the \"filter\" job successfully!");
			}		
		}
		else if(args[0].equals("barcode")){
			n = 3;
			String barCode = null;
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else if(args[i].equals("-b")){
					  barCode = args[i + 1];
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
			  if(barCode == null){
                  help();
				  System.exit(1);
			  }
		      SequenceUtils.setFormatLen(formatLen);
			  runBarCode(inputFileName, outputFileName, barCode);
			  System.out.println("Nector finished the \"barcode\" job successfully!");
			}

		}
	    else if(args[0].equals("exclude")){
			n = 3;
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
			  runExcluding(inputFileName, listFile, outputFileName);
			  System.out.println("Nector finished the \"exclude\" job successfully!");
			}
	    }
		else if(args[0].equals("select")){
			n = 3;
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
		      runSelect(inputFileName, listFile, outputFileName);
		      System.out.println("Nector finished the \"select\" job successfully!");
			}
		}
		else if(args[0].equals("fastq2fasta")){
			n = 3;
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
		      runConversion(inputFileName, outputFileName);
		      System.out.println("Nector finished the \"fastq2fasta\" job successfully!");
			}
		}
		else if(args[0].equals("sort")){
			n = 3;
			String order = new String("d");
            if(args.length < n){
               help();
			   System.exit(1);
			}
			else{
			  for(int i = n; i < args.length; i += 2){
                  if( (i + 1) == args.length){
                     help();
					 System.exit(1);
				  }
				  if(args[i].equals("-l")){
                     formatLen = Integer.parseInt(args[i + 1]);
				  }
				  else if(args[i].equals("-o")){
                     order = args[i + 1];
				  }
				  else{
                     help();
					 System.exit(1);
				  }
			  }//end for
		      SequenceUtils.setFormatLen(formatLen);
              runSortByLength(inputFileName, outputFileName, order);
              System.out.println("Nector finished the \"sort\" job successfully!");
			}
		}
		else{
             help();
		}
	}
	
    /**
	 * Converting FASTQ file to FASTA one
	 * @param inputFileName is the name of input FASTQ file
	 * @param outputFileName is the name of output FASTA file
	 */
	public static void runConversion(String inputFileName, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	      System.err.println("Error! Input file must be FASTQ format!");
		  System.exit(1);
		}
		else if(format == Format.FASTQ){
	      list = new SequenceListFastq(inputFileName);
          ((SequenceListFastq) list).convertToFasta(outputFileName);
        }
	}
	
    /**
	 * Trimming low quality Illumina reads and filtering short trimmed reads
	 * @param inputFileName is the name of input FASTQ file
	 * @param outputFileName is the name of output FASTA file
	 * @param cutoff is the minimum quality score, default value = 13
	 * @param maxPer is the maximum percentage of mismatches, defaul value = 0.05
	 * @param minLen is the minium trimmed length, default value = 40
	 */
	public static void runTrim(String inputFileName, String outputFileName, int cutoff, float  maxPer, int minLen){
        SequenceListFasta list = new SequenceListFastq(inputFileName);
        ((SequenceListFastq)list).trim(outputFileName, cutoff, maxPer, minLen);
	}
	
    /**
	 * Assembling the pair-ended reads into linked ones
	 * @param inputFileName is the name of input FASTQ file (forward)
	 * @param inputFileName2 is the name of input FASTQ file (reverse)
	 * @param outputFileName is the name of output FASTA file
	 * @param maxMismatch is the maximimum number of mismatches in overlapped region, default value = 3 
	 * @param minOverlap is the minimum bps of overlapped regions, default value = 10
	 */
	public static void runAssembly(String inputFileName, String inputFileName2, String outputFileName, int maxMismatch, int minOverlap){
		SequenceListFastaPair list = null;
		if(format == Format.FASTA){
           list = new SequenceListFastaPair(inputFileName + "," + inputFileName2);
		}
		else if(format == Format.FASTQ){
           list = new SequenceListFastqPair(inputFileName + "," + inputFileName2);
		}
		list.assembly(outputFileName, maxMismatch, minOverlap);
	}
	
    /**
	 * Giving statistics of FASTA/FASTQ file
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param outputFileName is the name of statistics output file
	 */
	public static void runStatistics(String inputFileName, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	      list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	      list = new SequenceListFastq(inputFileName);
		}
	    list.printSequenceStatistics(outputFileName);
	}
	
    /**
	 * Checking if the input file is FASTA or FASTQ file
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 */
	public static void runCheck(String inputFileName){
		SequenceListFasta list = null;
		boolean flag = true;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
           flag = list.formatCheck();
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
           flag = ((SequenceListFastq) list).formatCheck();
		}
		if(flag){
           if(format == Format.FASTA){
              System.out.println(inputFileName + " is a FASTA file.");
		   }
		   else if(format == Format.FASTQ){
              System.out.println(inputFileName + " is a FASTQ file.");
		   }
		}
		else{
           System.err.println("This file is NOT a FASTA or FASTQ file!");
		}
	}
	
    /**
	 * Get the sequences in reversed order
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param outputFileName is the name of output FASTA/FASTQ file in reverse order
	 */
	public static void runReverseOrder(String inputFileName, String outputFileName){
	    SequenceListFasta list = null;
		if(format == Format.FASTA){
		   list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
		   list = new SequenceListFastq(inputFileName);
		}
	    list.reverseSequence(outputFileName);
	}
	
    /**
	 * Get the sequences in complementary reverse
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param outputFileName is the name of output FASTA/FASTQ file in complementary reverse
	 */
	public static void runComplementaryReverse(String inputFileName, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
	    list.complementaryReverseSequence(outputFileName);
    }
	
    /**
	 * Print the sequences by specified line length for sequences
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param outputFileName is the name of output FASTA/FASTQ file by specified line length for sequences
	 */
	public static void runPrint(String inputFileName, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.printAllSequenceWithCache(outputFileName);
	}

    /**
	 * Get the sequences by ascending or descending order of length
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param outputFileName is the name of output FASTA/FASTQ file by ascending or descending order of length
     * @param order by ascending('a') or descending('d') order
	 */
	public static void runSortByLength(String inputFileName, String outputFileName, String order){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.sortSequenceByLen(outputFileName, order);
	}
	
    /**
	 * Filtering reads with shorter length
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param outputFileName is the name of output FASTA/FASTQ file (filtered reads)
     * @param minLen is the minimum length after filtering
	 */
	public static void runFilter(String inputFileName, String outputFileName, int minLen){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.filterSequence(outputFileName, minLen);
	}
	
    /**
	 * Filtering reads in a list file
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param listFile is the name of list containing the sequence description for excluding
	 * @param outputFileName is the name of output FASTA/FASTQ file (filtered reads)
	 */
	public static void runExcluding(String inputFileName, String listFile, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.excluding(listFile, outputFileName, true);
	}
	
    /**
	 * Extracting reads matching the barcode (barcode not included in the extracted sequences)
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param outputFileName is the name of output FASTA/FASTQ file (extracted reads)
     * @param barcode is the barcode information
	 */
	public static void runBarCode(String inputFileName, String outputFileName, String barcode){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.barCode(outputFileName, barcode);
	}
	
    /**
	 * Selecting reads in a list file
	 * @param inputFileName is the name of input FASTA/FASTQ file
	 * @param listFile is the name of list containing the sequence description for selected reads
	 * @param outputFileName is the name of output FASTA/FASTQ file (selected reads)
	 */
	public static void runSelect(String inputFileName, String listFile, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.excluding(listFile, outputFileName, false);
	}
	
    /**
	 * Fast check for input files 
	 * @param fileName is the name of input FASTA/FASTQ file
	 * @return the format of input file: FASTA/FASTQ/OTHER
	 */
	public static Format fastCheck(String fileName){
		String line;
		Format format = null;
		try{
	      Scanner in = new Scanner(new BufferedReader(new FileReader(fileName)));	
          if(in.hasNextLine()){
		     line = in.nextLine();
             if(line.substring(0, 1).equals(">")){
                format = Format.FASTA; 
			 }
             else if(line.substring(0, 1).equals("@")){
                format = Format.FASTQ; 
			 }
			 else{
                format = Format.OTHER;
			 }
          }
		  else{
                format = Format.OTHER;
		  }
		  in.close();
		}
		catch(IOException e){
		   e.printStackTrace();
	    }
		return format;
    }
}
