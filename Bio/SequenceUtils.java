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
 * <code>SequenceUtils</code> is a class for dealing with
 * static variables like sequence lenth of a line in FASTA
 * and FASTQ format.
 * @author Sitao Wu
 * @version 0.1
 * @since  2013
 */

public class SequenceUtils{
    private static int formatLen = 250;//default value = 250
    
	/**
     * Set sequence length in a line
     * @param len sequence length for a line 
	 */
	public static void setFormatLen(int len){
       formatLen = len;
	}
	
	/**
     * Get sequence length in a line
     * @return len sequence length for a line 
	 */
	public static int getFormatLen(){
       return formatLen;
	}

	/**
     * Get sequence length in a line
	 % @param  input a sequence or quality score line
     * @return multiple lines according to the formatLen for a line
	 */
	public static String format(String input){
	   String str = "";
	   int len = input.length();
       int m = len / formatLen;
	   int n = len % formatLen;
	   for(int i = 0; i < m; ++i){
	     str += input.substring(i * formatLen, (i + 1) * formatLen) + "\n"; 
	   }
	   if(n != 0) str += input.substring(m * formatLen, len) + "\n";
       return str; 
	}
	
	/**
	 * Asssembly two string together if there is overlapping
     * @param str1 the first string to be assembled
     * @param str2 the second string to be assembled
	 * @return assembly strings
	 */
	public static String [] doAssembly(String str1, String str2, int maxMismatch, int minOverlap){
		 int len1= str1.length();
		 int len2= str2.length();
		 int len = len1>len2?len2:len1;
		 int k = 0, tmp, i;
		 String desc = "";
		 char c1, c2;
		 boolean flag = false;
		 for(i = len; i > minOverlap ; --i){
		     flag = false;
			 k = 0;
			 desc = "";
			 for(int j = 0; j < i; ++j){
				tmp = len1 - i + j;
                c1 = str1.charAt(tmp);
                c2 = str2.charAt(j);
				if(c1 != c2){
                    ++k;
					if(k > maxMismatch){
					   flag = true;
                       break;
					}
					++tmp;
                    desc += " " + tmp + ":" + String.valueOf(c1) + ">" + String.valueOf(c2);
				}
			 }
			 if(!flag) break;
		 }
		 if(!flag){
			desc += " overlap_region=" + i + "bps"; 
		    return new String [] {desc, str1 + str2.substring(i, len2)};
		 }
		 else{
            return null;
		 }
	}
}
