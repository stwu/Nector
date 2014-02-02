Nector
======

Sequence preprocessing tool on NGS data


Program: Nector - sequence preprocessing tools for NGS data
Version 0.1. Developed by Sitao Wu
The MIT license

Usage:	java -jar Nector.jar toolName inputFileName(s) outputFileName parameters

Tools supported:

1.stat: Get statistics of sequence file

	Example 1: java -jar Nector.jar stat example.fasta stat.output
	Example 2: java -jar Nector.jar stat example.fastq stat.output

2.rev: Get sequences in revserse order

	Parameter: -l line width for sequence (default=250)

	Example 1: java -jar Nector.jar rev example.fasta exampleReverse.fasta [-l 50]
	Example 2: java -jar Nector.jar rev example.fastq exampleReverse.fastq [-l 50]

3.cr: Get sequences in complementary reverse

	Parameter: -l line width for sequence (default=250)

	Example 1: java -jar Nector.jar cr example.fasta exampleCr.fasta [-l 50]
	Example 2: java -jar Nector.jar cr example.fastq exampleCr.fastq [-l 50]

4.sort: Get sequences in sorted order (ascending or descending) by length

	Parameter: -l line width for sequence (default=250)
	Parameter: -o a:ascending;d:descending (default=d)

	Example 1: java -jar Nector.jar sort example.fasta exampleSort.fasta [-l 50] (default sort order = descending)
	Example 2: java -jar Nector.jar sort example.fastq exampleSort.fastq [-l 50] (default sort order = descending)
	Example 3: java -jar Nector.jar sort example.fasta exampleSort.fasta [-o a] [-l 50]
	Example 4: java -jar Nector.jar sort example.fastq exampleSort.fastq [-o a] [-l 50]
	Example 5: java -jar Nector.jar sort example.fasta exampleSort.fasta [-o d] [-l 50]
	Example 6: java -jar Nector.jar sort example.fastq exampleSort.fastq [-o d] [-l 50]

5.filter: filter short sequences

	Parameter: -l line width for sequence (default=250)
	Parameter: -n minimum length (default=40)

	Example 1: java -jar Nector.jar filter example.fasta exampleFiltered.fasta [-n 100] [-l 50] (minimum sequence length=100)
	Example 2: java -jar Nector.jar filter example.fastq exampleFiltered.fasta [-n 100] [-l 50] (minimum sequence length=100)

6.exclude: excluding sequences in a list

	Parameter: -l line width for sequence (default=250)

	Example 1: java -jar Nector.jar exclude example.fasta,list exampleExcluded.fasta [-l 50] (excluding reads defined in list by read name)
	Example 2: java -jar Nector.jar exclude example.fastq,list exampleExcluded.fastq [-l 50] (excluding reads defined in list by read name)

7.select: selecting sequences in a list

	Parameter: -l line width for sequence (default=250)

	Example 1: java -jar Nector.jar select example.fasta,list exampleSelecteded.fasta [-l 50] (selecting reads defined in list by read name)
	Example 2: java -jar Nector.jar select example.fastq,list exampleSelecteded.fastq [-l 50] (selecting reads defined in list by read name)

8.fastq2fasta: converting fastq to fasta format

	Parameter: -l line width for sequence (default=250)

	Example: java -jar Nector.jar fastq2fasta example.fastq exampleConvert.fasta [-l 50]

9.check: check format

	Example 1: java -jar Nector.jar check example.fasta
	Example 2: java -jar Nector.jar check example.fastq

10.trim: trimming and filtering

	Parameter: -l line width for sequence (default=250)
	Parameter: -c minimum quality score (default=13)
	Parameter: -p maximum error rate (default=0.05)
	Parameter: -L minimum sequence length (default=40

	Example: java -jar Nector.jar trim example.fastq exampleTrimed.fastq [-c 13] [-p 0.05] [-L 70] [-l 50]

11.print: write file according to specified lengh per line

	Parameter: -l line width for sequence (default=250)

	Example: java -jar Nector.jar print example.fasta exampleNew.fasta [-l 50]
	Example: java -jar Nector.jar print example.fastq exampleNew.fastq [-l 50]

12.barcode: extract sequences according to barcode

	Parameter: -l line width for sequence (default=250)
	Parameter: -b bar code to be filtered

	Example: java -jar Nector.jar barcode example.fasta exampleBarcode.fasta -b ACC [-l 50]
	Example: java -jar Nector.jar barcode example.fastq exampleBarcode.fastq -b ACC [-l 50]

13.assembly: assembly pair-ended reads into linked ones

	Parameter: -l line width for sequence (default=250)
	Parameter: -m maximum number of mismatch (default=3)
	Parameter: -v minimum overlapping length (default=10)

	Example: java -jar Nector.jar assembly example1.fasta,example2.fasta exampleAssembled.fasta [-m 3] [-v 10] [-l 50]
	Example: java -jar Nector.jar assembly example1.fastq,example2.fastq exampleAssembled.fasta [-m 3] [-v 10] [-l 50]

