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


/**
 * An object of <code>DnaSequenceFasta</code> is a 
 * dna sequence file in FASTA format
 * @author Sitao Wu
 * @version 0.1
 * @since  2013
 */
public class DnaSequenceFasta extends SequenceFasta {
	private static char [] complement = { 'T', 'B', 'G', 'D', 'E', 'F', 'C',
	                                      'H', 'I', 'J', 'K', 'L', 'M', 'N',
							              'O', 'P', 'Q',      'R', 'S', 'A',
								          'A', 'V', 'W',      'X', 'Y', 'Z'};


    /**
     * DnaSequenceFasta constructor
	 */
    public DnaSequenceFasta(){
        super();
	}

    /**
     * DnaSequenceFasta constructor
	 * @param description for a desciption line
	 * @param sequence for a sequence line 
	 */
    public DnaSequenceFasta(String description, String sequence){
        super(description, sequence);
	}

    /**
     * DnaSequenceFasta constructor
	 * @param fasta for a DnaSequenceFasta object
	 */
    public DnaSequenceFasta(DnaSequenceFasta fasta){
        super(fasta); 
	}

    /**
     * Check if the sequence is a DNA sequence
	 * @return true if the sequence is a DNA sequence; false otherwise
	 */
	public boolean sequenceCheck(){
	   char [] seq = getSequence().toCharArray();
	   int [] count = new int[26];
	   int i, len;
	   for(i = 0; i < 26; i++){
          count[i] = 0;
	   }
	   len = seq.length;
       for(char a: seq){
		  count[a - 'A']++; 
	   }
       int sum = 0;
	   char [] tmp = {'A', 'C', 'G', 'T', 'N'};
	   for(char t: tmp){
		   sum += count[t - 'A'];   
	   }
	   if(sum == len){
	      return true;
	   }
	   else{
		  System.err.println("Error! Sequences should only contains A, C, G, T, N");
          return false;
	   }
	}
 
    /**
     * Check if the description begins with a ">" symbol
	 * @return true if the description begins with a ">" symbol; false otherwise
	 */
	public boolean descriptionCheck(){
       if(!getDescription().substring(0,1).equals(">")){
		    System.err.println("Error! Description should begin with \">\" in FASTA format!");
            return false;
	   }
	   else{
            return true;
	   }
	}
	
    /**
     * Check if the DnaSequenceFasta is a FASTA format using sequenceCheck() and descriptionCheck() 
	 * @return true if the DnaSequenceFasta is a FASTA format; false otherwise
	 */
	public boolean formatCheck(){
	   //check quality area
       if(!descriptionCheck()){
            return false;
	   }
	   //check sequence
       if(!sequenceCheck()){
            return false;
	   }
	   return true;
	}
	
    /**
     * Make a complementary reverse sequence from the original DnaSequenceFasta object 
	 * @return a complementary reverse sequence in FASTA format
	 */
	public DnaSequenceFasta complementaryReverse(){
        int i, j, len, mid;
		int left, right;
		char [] array;
		len = getLength();
		array = getSequence().toCharArray();
		mid = (len / 2) + (len % 2);
		for(i = 0, j = len - 1; i < mid; i++, j--){
           left = array[i];
		   right = array[j];
           array[i] = complement[right - 'A']; 
           array[j] = complement[left - 'A']; 
		}
		return new DnaSequenceFasta(getDescription(), new String(array));
	}
	
	/**
	 * Use barcode to extract sequences
     * @param barCode barcode sequence
	 * @return extracted sequences. Null if barcode does not match
	 */
	 public DnaSequenceFasta barcode(String barCode){
        int len = barCode.length();
		String str = getSequence();
		int len2 = str.length();
		if(str.substring(0, len).equals(barCode)){
           return new DnaSequenceFasta(getDescription(), str.substring(len, len2));      
		}
		else{
           return null;
		}
	 }
	 
	/**
	 * Assembly two overlapped pair-ended reads together into current read
     * @param other the 3' end read data
     * @param maxMismatch the maximum number of mismatches in overlapped region
     * @param minOverlap the minimum length of overlapped region
     * @return true if assembly is successful; false otherwise
	 */
	 public boolean assembly(DnaSequenceFasta other, int maxMismatch, int minOverlap){
        DnaSequenceFasta cr = other.complementaryReverse();
		String [] str = SequenceUtils.doAssembly(getSequence(), cr.getSequence(), maxMismatch, minOverlap);
		if(str != null){
		   setDescription(getDescription() + " " + str[0]);
		   setSequence(str[1]);
		   return true;
		}
		else{
           return false;
		}
	 }	 
}
