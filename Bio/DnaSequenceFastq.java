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
 * An object of <code>DnaSequenceFastq</code> is a 
 * dna sequence file in FASTQ format
 * @author Sitao Wu
 * @version 0.1
 * @since  2013
 */
public class DnaSequenceFastq extends DnaSequenceFasta{
   private enum Platform{
	   NONE,
	   SANGER,
	   ILLUMINA_1_0,
	   ILLUMINA_1_3 
   }
   private String qDescription;//Quality description
   private String quality;//Quality score
   private Platform platform;//platform

   /**
    * DnaSequenceFastq constructor
    */
   public DnaSequenceFastq(){
	  super();
	  qDescription = null;
	  quality = null;
	  platform = Platform.ILLUMINA_1_3; 
   }
   
   /**
    * DnaSequenceFastq constructor
	* @param fastq is an object of DnaSequenceFastq
    */
   public DnaSequenceFastq(DnaSequenceFastq fastq){
      setFastq(fastq);
   }
   
   /**
    * DnaSequenceFastq constructor
	* @param description is a description line
	* @param sequence is a sequence line
	* @param qDescription is the description line for quality sequence
	* @param quality is a quality score line
	* @param platform is the platform description
    */
   public DnaSequenceFastq(String description, String sequence, String qDescription, String quality, Platform platform){
      this(description, sequence, qDescription, quality);
      setPlatform(platform);
   }
   /**
    * DnaSequenceFastq constructor
    * DnaSequenceFastq constructor
	* @param description is a description line
	* @param sequence is a sequence line
	* @param qDescription is the description line for quality sequence
	* @param quality is a quality score line
    */
   public DnaSequenceFastq(String description, String sequence, String qDescription, String quality){
	  super(description, sequence);
	  this.qDescription = qDescription;
	  this.quality = quality;
      setPlatform(platformCheck());
   }
   
   /**
     * Set the quality description line
	 * @param qDescription is the quality description line
	 */
   public void setQdescription(String qDescription){
	  this.qDescription = qDescription;
   }
   
   /**
     * Set the quality score line
	 * @param quality is the quality score line
	 */
   public void setQuality(String quality){
	  this.quality = quality;
   }
   
   /**
     * Set the platform
	 * @param platform is the platform description
	 */
   public void setPlatform(Platform platform){
	  this.platform = platform;
   }
   
   /**
     * Get the quality description line 
	 * @return the quality description line
	 */
   public String getQdescription(){
      return qDescription;
   }
   
   /**
     * Get the quality score line 
	 * @return the quality score line
	 */
   public String getQuality(){
      return quality;
   }
   
   /**
     * Get the platform description 
	 * @return the platform description
	 */
   public Platform getPlatform(){
      return platform;
   }
   
   /**
     * Trim and filter the sequence according to quality score
	 * @param cutoff the minimum quality score. Quality score=13 (default value) is corresponding to error rate= 5%. 
	 * @param maxPer the maximum bases below cutoff. Default value = 0.05
	 * @param minLen the minimum len after trimming. If length after trimming < minLen, set length to zero. Default value = 40
	 */
   public void trim(int cutoff, float maxPer, int minLen){
       String str;
	   str = getSequence();
       char [] qual = quality.toCharArray();
	   int i = 0, j = 0;
	   int len = quality.length();
       for(char s: qual){
		 j++;
         if((s - '@') < cutoff){
            i++;
            float tmp = (float) i / j;
			if(tmp >= maxPer){
				if(i > minLen){
                   setSequence(str.substring(0, i - 1));
				   setQuality(quality.substring(0, i - 1));
				}
				else{//remove this sequence
                   setSequence("");
                   setQuality("");
				}
				return;
			}
		 }
	   }
   }
   
   /**
     * Check if the description begins with a "@" symbol
	 * @return true if the description begins with a "@" symbol; false otherwise
	 */
   @Override
   public boolean descriptionCheck(){
       if(!getDescription().substring(0,1).equals("@")){
		    System.err.println("Error! Sequence description should begin with \"@\" in FASTQ format!");
            return false;
	   }
	   else{
            return true;
	   }
   }
   
   /**
     * Check if the quality description begins with a "+" symbol
	 * @return true if the description begins with a "+" symbol; false otherwise
	 */
   public boolean qDescriptionCheck(){
       if(!qDescription.substring(0,1).equals("+")){
		    System.err.println("Error! Quality description should begin with \"+\" in FASTQ format!");
            return false;
	   }
	   else{
            return true;
	   }
   }
   
   /**
     * Check the platform of the FASTQ file
	 * @return the platform of the FASTQ file
	 */
   public Platform platformCheck(){
      char [] seq = quality.toCharArray();
	  boolean bNone, bSanger, bIllumina_1_0, bIllumina_1_3;
	  char high = '!', low = '~';
	  
      for(char s : seq){
         if( s < low) low = s;
		 if( s > high) high = s; 		  
	  }
	  if( low >= '@' && high <= '~'){
         platform = Platform.ILLUMINA_1_3;
	  }
	  else if( low >= ';' && high <= '~'){
         platform = Platform.ILLUMINA_1_0;
	  }
	  else if( low >= '!' && high <= '~'){
         platform = Platform.SANGER;
	  }
	  else{
         platform = Platform.NONE;
		 System.err.println("Error! Quality score out of range [33-126]!");
	  }
	  return platform;
   }
   
   /**
     * Check sequence and quality have the same length
	 * @return true if sequence and quality have the same length; false otherwise
	 */
   public boolean lengthCheck(){
       if(getLength() != quality.length()){//sequence length is the same as that of quality scores
		  System.err.println("Error! The lengths of sequence and quality score are different!");
          return false;
	   }
	   else{
		   return true;
	   }
   }
   
   /**
     * Check the FASTQ format using descriptionCheck(), sequenceCheck(), qdescriptionCheck(), platformCheck() and lengthCheck() 
	 * @return true if FASTQ format is correct; false otherwise
	 */
   public boolean formatCheck(){
       if(!descriptionCheck()){//description  check
          return false;
	   }
       if(!sequenceCheck()){//sequence  check
          return false;
	   }
       if(!qDescriptionCheck()){//Qdescription check
	      return false;
	   }
	   if(platformCheck() == Platform.NONE){//platform check
          return false;
	   }
       if(!lengthCheck()){//sequence length is the same as that of quality scores
          return false;
	   }
       return true;
   }
   
   /**
    * Set a DnaSequenceFastq object
    * @param fastq a DnaSequenceFastq object
    */
   public void setFastq(DnaSequenceFastq fastq){
         setFasta(fastq.convertToFasta(fastq));
		 qDescription = fastq.getQdescription();
		 quality = fastq.getQuality();
		 platform = fastq.getPlatform();
   }

   /**
     * Give the print content including description, sequence, quality description, and quality score
     * @return a String object containing the description, sequence, quality description, and quality score
	 */
   @Override
   public String toString(){
       return getDescription() + "\n" + SequenceUtils.format(getSequence()) + qDescription + "\n" + SequenceUtils.format(quality);
   }

   /**
    * Convert fastq to fasta format
    * @return a DnaSequenceFasta object
    */
   public DnaSequenceFasta convertToFasta(){
	 String tmpStr = getDescription();
	 tmpStr = ">" + tmpStr.substring(1, tmpStr.length());
	 return new DnaSequenceFasta(tmpStr, getSequence());  
   }
   
   /**
    * Convert fastq to fasta format. Static function
    * @param fastq a DnaSequenceFastq object
    * @return a DnaSequenceFasta object
    */
   public static SequenceFasta convertToFasta(DnaSequenceFastq fastq){
	 return fastq.convertToFasta();  
   }
   
   /**
    * Make a complementary reverse sequence from the original DnaSequenceFastq object 
    * @return a complementary reverse sequence in FASTQ format
	*/
   @Override
   public DnaSequenceFasta complementaryReverse(){
     DnaSequenceFasta fasta = convertToFasta();
	 fasta = fasta.complementaryReverse();
	 int i, j;
	 char tmp;
	 int len = getLength();
	 char [] array = getQuality().toCharArray();
 	 int mid = (len / 2) + (len % 2);
	 for(i = 0, j = len - 1; i < mid; i++, j--){
           tmp = array[i];
           array[i] = array[j]; 
           array[j] = tmp;
	 }
	 String tmpStr = getDescription();
	 tmpStr = "@" + tmpStr.substring(1, tmpStr.length());
	 return new DnaSequenceFastq(tmpStr, fasta.getSequence(), getQdescription(), new String(array));
   }
   
   /**
	 * Use barcode to extract sequences
     * @param barCode barcode sequence
	 * @return extracted sequences if successful; null otherwise
	 */
   public DnaSequenceFastq barcode(String barCode){
     int len = barCode.length();
	 String str = getSequence();
	 int len2 = str.length();
	 if(str.substring(0, len).equals(barCode)){
        return new DnaSequenceFastq(getDescription(), str.substring(len, len2), getQdescription(), getQuality().substring(len, len2));      
	 }
	 else{
        return null;
	 }
   }
}
