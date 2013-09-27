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
        System.out.println("Usage:\tNector toolName fileName [parameters]");
		System.out.println("");
		System.out.println("Tools supported:");
		System.out.println("");
		System.out.println("1.stat: print statistics of sequence file");
		System.out.println("");
		System.out.println("\tExample: Nector stat test.fasta");
		System.out.println("");
		System.out.println("2.rc: print sequences in revserse complementary");
		System.out.println("");
		System.out.println("\tExample: Nector rc test.fasta");
		System.out.println("");
	}
	public static void processArg(String [] args){
		String fileName;
        if(args.length == 0){
            help();
			return;
		}
		if(args[0].equals("stat")){
            if(args.length != 2){
			   help();
			   return;	
			}
			fileName = args[1];
			runStatistics(fileName);
		}
		else if(args[0].equals("rc")){
            if(args.length != 2){
			   help();
			   return;	
			}
			fileName = args[1];
			runReverseComplement(fileName);
		}
	}
	public static void runStatistics(String fileName){
	    SequenceList list = new SequenceList(fileName);
	    list.printSequenceStatistics();
	    //list.sortSequenceByLen();
	}
	public static void runReverseComplement(String fileName){
	    SequenceList list = new SequenceList(fileName);
	    list.reverseSequence();
	    list.printAllSequences();
        
	}
}
