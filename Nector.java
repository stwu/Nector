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

public class Nector{
	 private enum Format{
		 OTHER,
		 FASTA,
		 FASTQ
    };
    private static Format format;
	private static String version = "0.1";
    public static void main(String [] args){
		processArg(args);
	}
	public static void help(){
		System.out.println("");
		System.out.println("Program: Nector - sequence processing tools");
		System.out.printf("Version %s. Developed by Sitao Wu\n", version);
		System.out.println("");
        System.out.println("Usage:\tjava Nector toolName inputFileName outputFileName [parameters]");
		System.out.println("");
		System.out.println("Tools supported:");
		System.out.println("");
		System.out.println("1.stat: Get statistics of sequence file");
		System.out.println("");
		System.out.println("\tExample 1: java Nector stat example.fasta stat.txt");
		System.out.println("\tExample 2: java Nector stat example.fastq stat.txt");
		System.out.println("");
		System.out.println("2.rev: Get sequences in revserse order");
		System.out.println("");
		System.out.println("\tExample 1: java Nector rev example.fasta exampleReverse.fasta");
		System.out.println("\tExample 2: java Nector rev example.fastq exampleReverse.fastq");
		System.out.println("");
		System.out.println("3.cr: Get sequences in complementary reverse");
		System.out.println("");
		System.out.println("\tExample 1: java Nector cr example.fasta exampleCr.fasta");
		System.out.println("\tExample 2: java Nector cr example.fastq exampleCr.fastq");
		System.out.println("");
		System.out.println("4.sort: Get sequences in sorted order (ascending or descending) by length");
		System.out.println("");
		System.out.println("\tExample 1: java Nector sort example.fasta exampleSort.fasta (default sort order = descending)");
		System.out.println("\tExample 2: java Nector sort example.fastq exampleSort.fastq (default sort order = descending)");
		System.out.println("\tExample 3: java Nector sort example.fasta exampleSort.fasta ascending");
		System.out.println("\tExample 4: java Nector sort example.fastq exampleSort.fastq ascending");
		System.out.println("\tExample 5: java Nector sort example.fasta exampleSort.fasta descending");
		System.out.println("\tExample 6: java Nector sort example.fastq exampleSort.fastq descending");
		System.out.println("");
		System.out.println("5.filter: filter short sequences");
		System.out.println("");
		System.out.println("\tExample 1: java Nector filter example.fasta exampleFiltered.fasta 100 (minimum sequence length=100)");
		System.out.println("\tExample 2: java Nector filter example.fastq exampleFiltered.fasta 100 (minimum sequence length=100)");
		System.out.println("");
		System.out.println("6.exclude: excluding sequences in a list");
		System.out.println("");
		System.out.println("\tExample 1: java Nector exclude example.fasta,list exampleExcluded.fasta (excluding reads defined in list by read name)");
		System.out.println("\tExample 2: java Nector exclude example.fastq,list exampleExcluded.fasta (excluding reads defined in list by read name)");
		System.out.println("");
		System.out.println("7.select: selecting sequences in a list");
		System.out.println("");
		System.out.println("\tExample 1: java Nector select example.fasta,list exampleSelecteded.fasta (selecting reads defined in list by read name)");
		System.out.println("\tExample 2: java Nector select example.fastq,list exampleSelecteded.fasta (selecting reads defined in list by read name)");
		System.out.println("");
		System.out.println("8.fastq2fasta: converting fastq to fasta format");
		System.out.println("");
		System.out.println("\tExample: java Nector fastq2fasta example.fastq exampleConvert.fasta");
		System.out.println("");
		System.out.println("9.check: check format");
		System.out.println("");
		System.out.println("\tExample 1: java Nector check example.fasta");
		System.out.println("\tExample 2: java Nector check example.fastq");
		System.out.println("");
		System.out.println("10.trim: trimming and filtering");
		System.out.println("");
		System.out.println("\tExample: java Nector trim example.fastq exampleTrimed.fastq 13 0.05 70");
		System.out.println("");
	}
	public static void processArg(String [] args){
		String inputFileName = null, listFile = null;
		String outputFileName = null;
        if(args.length < 2){
		  help();
		  return;
		}
        else if(args.length == 2 && !(args[0].equals("check"))){
		  help();
		  return;	
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
		format = fastCheck(inputFileName);
		if(format == Format.OTHER){
			 System.err.println("Error! Input file is not FASTA or FASTQ file!");
             System.exit(1);
		}
		if(args[0].equals("stat")){
			runStatistics(inputFileName, outputFileName);
		}
		else if(args[0].equals("check")){
            runCheck(inputFileName);
		}
		else if(args[0].equals("rev")){
			runReverseOrder(inputFileName, outputFileName);
		}
		else if(args[0].equals("cr")){
			runComplementaryReverse(inputFileName, outputFileName);
		}
		else if(args[0].equals("trim")){
            if(args.length < 6){
               help();
			}
			else{
               int cutoff = Integer.parseInt(args[3]);
			   float maxPer = Float.parseFloat(args[4]);
			   int minLen = Integer.parseInt(args[5]);
			   runTrim(inputFileName, outputFileName, cutoff, maxPer, minLen);
			}
		}
		else if(args[0].equals("filter")){
            if(args.length == 3){
               help();
			}
            else if(args.length == 4){
			   int minLen = 0;
			   try{
                   minLen = Integer.parseInt(args[3]);
			   }
			   catch(NumberFormatException e){
				   System.err.println("Error! Please give an integer for a minimum length cutoff");
				   System.exit(1);
			   }
			   runFilter(inputFileName, outputFileName, minLen);
			}
		}
	    else if(args[0].equals("exclude")){
			runExcluding(inputFileName, listFile, outputFileName);
	    }
		else if(args[0].equals("select")){
		   runSelect(inputFileName, listFile, outputFileName);
		}
		else if(args[0].equals("fastq2fasta")){
		   runConversion(inputFileName, outputFileName);
		}
		else if(args[0].equals("sort")){
            if(args.length == 4){
               if(args[3].equals("ascending") || args[3].equals("descending")){
                  runSortByLength(inputFileName, outputFileName, args[3]);
			   }
			   else{
                  help();
			   }
			}
			else{
               runSortByLength(inputFileName, outputFileName, "descending");
			}
		}
		else{
             help();
		}
	}
	public static void runConversion(String inputFileName, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	      System.err.println("Error! Input file must be FASTQ format!");
		  System.exit(1);
		}
		else if(format == Format.FASTQ){
	      list = new SequenceListFastq(inputFileName);
          ((SequenceListFastq) list).convertToFasta();
		  ((SequenceListFastq) list).printAllSequencesFasta(outputFileName);
        }
	}
	public static void runTrim(String inputFileName, String outputFileName, int cutoff, float  maxPer, int minLen){
        SequenceListFasta list = new SequenceListFastq(inputFileName);
        ((SequenceListFastq)list).trim(cutoff, maxPer, minLen);
	    list.printAllSequences(outputFileName);	 
	}
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
              System.out.println(inputFileName + " is a FASTA file");
		   }
		   else if(format == Format.FASTQ){
              System.out.println(inputFileName + " is a FASTQ file");
		   }
		}
		else{
           System.err.println("This file is NOT a FASTA or FASTQ file!");
		}

	}
	public static void runReverseOrder(String inputFileName, String outputFileName){
	    SequenceListFasta list = null;
		if(format == Format.FASTA){
		   list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
		   list = new SequenceListFastq(inputFileName);
		}
	    list.reverseSequence();
	    list.printAllSequences(outputFileName);
	}
	public static void runComplementaryReverse(String inputFileName, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
	    list.complementaryReverseSequence();
	    list.printAllSequences(outputFileName);
    }
	public static void runSortByLength(String inputFileName, String outputFileName, String order){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.sortSequenceByLen(order);
	    list.printAllSequences(outputFileName);
	}
	public static void runFilter(String inputFileName, String outputFileName, int minLen){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.filterSequence(minLen);
	    list.printAllSequences(outputFileName);
	}
	public static void runExcluding(String inputFileName, String listFile, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.excluding(listFile, true);
	    list.printAllSequences(outputFileName);
	}
	public static void runSelect(String inputFileName, String listFile, String outputFileName){
		SequenceListFasta list = null;
		if(format == Format.FASTA){
	       list = new SequenceListFasta(inputFileName);
		}
		else if(format == Format.FASTQ){
	       list = new SequenceListFastq(inputFileName);
		}
		list.excluding(listFile, false);
	    list.printAllSequences(outputFileName);
	}
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
