import Bio.*;
import java.io.*;
import java.util.*;

public class Nector{
	private static String version = "0.1";
    public static void main(String [] args){
		processArg(args);
	}
	public static void help(){
		System.out.printf("Nector: NGS data processing tools, version %s. Developed by Sitao Wu, email:wusitao2000@gmail.com\n", version);
		System.out.println("");
        System.out.println("Usage:\tjava Nector toolName inputFileName outputFileName [parameters]");
		System.out.println("");
		System.out.println("Tools supported:");
		System.out.println("");
		System.out.println("1.stat: Get statistics of sequence file");
		System.out.println("");
		System.out.println("\tExample: java Nector stat test.fasta stat.txt");
		System.out.println("");
		System.out.println("2.rev: Get sequences in revserse order");
		System.out.println("");
		System.out.println("\tExample: java Nector rev test.fasta testReverse.fasta");
		System.out.println("");
		System.out.println("3.cr: Get sequences in complementary reverse");
		System.out.println("");
		System.out.println("\tExample: java Nector cr test.fasta testCr.fasta");
		System.out.println("");
		System.out.println("4.sort: Get sequences in sorted order (ascending or descending) by length");
		System.out.println("");
		System.out.println("\tExample 1: java Nector sort test.fasta testSort.fasta (default sort order = descending)");
		System.out.println("\tExample 2: java Nector sort test.fasta testSort.fasta ascending");
		System.out.println("\tExample 3: java Nector sort test.fasta testSort.fasta descending");
		System.out.println("");
		System.out.println("5.filter: filter short sequences");
		System.out.println("");
		System.out.println("\tExample: java Nector filter te$st.fasta testFiltered.fasta 100 (minimum sequence length=100)");
		System.out.println("");
	}
	public static void processArg(String [] args){
		String inputFileName;
		String outputFileName;
        if(args.length < 3){
		  help();
		  return;	
		}
		inputFileName = args[1];
		outputFileName = args[2];

		if(args[0].equals("stat")){
			runStatistics(inputFileName, outputFileName);
		}
		else if(args[0].equals("rev")){
			runReverseOrder(inputFileName, outputFileName);
		}
		else if(args[0].equals("cr")){
			runComplementaryReverse(inputFileName, outputFileName);
		}
		else if(args[0].equals("filter")){
            if(args.length > 4){
               help();
			}
            if(args.length == 4){
			   int minLen = 0;
			   try{
                   minLen = Integer.parseInt(args[3]);
			   }
			   catch(NumberFormatException e){
				   System.out.println("Error! Please give an integer for a minimum length cutoff");
				   System.exit(-1);
			   }
			   runFilter(inputFileName, outputFileName, minLen);
			}
		}
		else if(args[0].equals("sort")){
            if(args.length > 4){
               help();
			}
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
	public static void runStatistics(String inputFileName, String outputFileName){
	    SequenceListFasta list = new SequenceListFasta(inputFileName);
	    list.printSequenceStatistics(outputFileName);
	    //list.sortSequenceByLen();
	}
	public static void runReverseOrder(String inputFileName, String outputFileName){
	    SequenceListFasta list = new SequenceListFasta(inputFileName);
	    //list.printAllSequences("before.txt");
	    list.reverseSequence();
	    list.printAllSequences(outputFileName);
	}
	public static void runComplementaryReverse(String inputFileName, String outputFileName){
	    SequenceListFasta list = new SequenceListFasta(inputFileName);
	    list.complementaryReverseSequence();
	    list.printAllSequences(outputFileName);
    }
	public static void runSortByLength(String inputFileName, String outputFileName, String order){
	    SequenceListFasta list = new SequenceListFasta(inputFileName);
		list.sortSequenceByLen(order);
	    list.printAllSequences(outputFileName);
	}
	public static void runFilter(String inputFileName, String outputFileName, int minLen){
	    SequenceListFasta list = new SequenceListFasta(inputFileName);
		list.filterSequence(minLen);
	    list.printAllSequences(outputFileName);
	}
}
