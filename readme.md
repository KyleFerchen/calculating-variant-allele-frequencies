# Example implementation for calculating variant allele frequencies

These scripts try to set up a method for calculating variant allele frequencies (VAFs) from mitochondrial reads of a 10X Genomics scRNA-seq dataset.

### Isolating mitochondrial reads

First, the output bam file from the 10X Genomics pipeline is filtered to only the reads that mapped to the mitochondrial genome.
You can check the name of the chromosome to filter by using the samtools view command.
Example:
`samtools view -H in.bam`

Then, use the name provided for the mitochondrial genome to filter to the reads mapped to that chromosome.
Example:
`samtools view -b in.bam chrMT > out.bam`

### Counting nucleotide occurrences

Next, the script 'get_nucleotide_counts_for_mt_reads.py' can be run on the filtered bam file to count the occurrences for each base at each position that has coverage from the experiment. The result is stored in a shelf file as a dictionary.

### Calculating allele frequencies

Lastly, a script can be run that takes the dictionary created from the previous step, and calculates the frequency for each allele.