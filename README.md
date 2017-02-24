# TelomereCounter
Counts telomeric repeats in the whole genome sequencing fastqs in GZIP textfiles 

Inputs: 
	 * args[0]: motif (e.g. TTAGGG)
	 * args[1]: min # of total telomeric repeats in paired end reads (e.g. 6 with at least 3 repeats per read)
	 * args[2]: FASTQ GZip textfile of telomeric reads 
	 * args[3]: directory for outputs

Example input 1: 
java -jar TelomereCounter.jar TTAGGG 6 /path/to/readFile.txt.gz /path/to/Output_directory/

Example input 2:
java -jar TelomereCounter.jar TTAGGG 6 /Users/kendrahong/Desktop/Lab/JAR/Test/Reads/SRR1008518_2.txt.gz /Users/kendrahong/Desktop/Lab/JAR/Test/Output/
