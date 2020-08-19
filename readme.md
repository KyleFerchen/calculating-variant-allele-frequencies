# Example implementation for calculating variant allele frequencies

These scripts try to set up a method for calculating variant allele frequencies (VAFs) from mitochondrial reads of a 10X Genomics scRNA-seq dataset. This is an attempt at following the method to perform lineage tracing using mitochondrial mutations and single-cell genomics described [here.](https://www.cell.com/cell/pdf/S0092-8674(19)30055-8.pdf)
<br>
[Ludwig, L. S., Lareau, C. A., Ulirsch, J. C., Christian, E., Muus, C., Li, L. H., ... & Law, T. (2019). Lineage tracing in humans enabled by mitochondrial mutations and single-cell genomics. Cell, 176(6), 1325-1339.](https://www.cell.com/cell/pdf/S0092-8674(19)30055-8.pdf)

### Isolating mitochondrial reads

First, the output bam file from the 10X Genomics pipeline is filtered to only the reads that mapped to the mitochondrial genome.
You can check the name of the chromosome to filter by using the samtools view command.
<br>
Example:
<br>
`samtools view -H in.bam`

Then, use the name provided for the mitochondrial genome to filter to the reads mapped to that chromosome.
<br>
Example:
<br>
`samtools view -b in.bam chrMT > out.bam`

### Counting nucleotide occurrences

Next, the script 'get_nucleotide_counts_for_mt_reads.py' can be run on the filtered bam file to count the occurrences for each base at each position that has coverage from the experiment. The result is stored in a shelf file as a dictionary. It uses pysam to loop through the bam reads, counting the occurrence for each base at each position.

### Calculating allele frequencies

Lastly, a script can be run that takes the dictionary created from the previous step, and calculates the frequency for each allele.

For example, 'count_af_h34.py' in this repository.