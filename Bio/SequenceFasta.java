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
 * An object of <code>SequenceFasta</code> is a 
 * general purpose sequence file in FASTA format
 * @author Sitao Wu
 * @version 0.1
 * @since  2013
 */
public class SequenceFasta implements Comparable<SequenceFasta> {
	private String description;//Desciption of sequences, treat it like ID, don't include ">"
    private String sequence;//Sequence like ACGT
    
    /**
     * SequenceFasta constructor
	 */
	public SequenceFasta(){
       description = null;
	   sequence = null;
	}

	/**
     * SequenceFasta constructor
	 * @param fasta for a SequenceFasta object
	 */
    public SequenceFasta(SequenceFasta fasta){
       setFasta(fasta);
	}
    
	/**
	 * SequenceFasta constructor
	 * @param description for a description line
	 * @param sequence for a sequence line
	 */
	public SequenceFasta(String description, String sequence){
	   this.description = description;
       this.sequence = sequence.toUpperCase();
	}
    
	/**
     * Set a sequence
	 * @param sequence for a sequence line
	 */
	public void setSequence(String sequence){
       this.sequence = sequence.toUpperCase();
	}
    
	/**
     * Set a description
     * @param description for a description line
	 */
	public void setDescription(String description){
       this.description = description;
	}
    
	/**
     * Set a SequenceFasa object
     * @param fasta for a SequenceFasta object
	 */
	public void setFasta(SequenceFasta fasta){
	   description = fasta.getDescription(); 
       sequence = fasta.getSequence();
	}

    /**
     * Get a sequence
     * @return a sequence line
	 */
	public String getSequence(){
       return sequence;
	}

    /**
     * Get a description
     * @return a description line
	 */
	public String getDescription(){
       return description;
	}
    
	/**
     * Get the length of sequence
     * @return the length of sequence
	 */
	public int getLength(){
       return sequence.length(); 
	}
    
	/**
     * Compare the length of this object to other one
     * @param other SequenceFasa object
     * @return negative or zero or positve if the length of this object < or = or > that of other object
	 */
	@Override
	public int compareTo(SequenceFasta other){
       return sequence.length() - other.sequence.length();
	}
    
	/**
     * Give the print content including description & sequence
     * @return a String object containing the description & sequence
	 */
	@Override
	public String toString(){
       return description + "\n" + SequenceUtils.format(sequence);
	}
}
