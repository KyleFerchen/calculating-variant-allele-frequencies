'''
This program takes in a bam file containing reads mapped to the mitochondrial genome from
10X genomics single cell RNA-Seq alignments (Contains 'CB' tag in reads). It returns a
python shelf file that stores a dictionary of cells (unique 'CB' tags), which each contain
a dictionary of nucleotides (A, T, C, G), each of which contains a dictionary of positions
and the number of times the base was read in that position. A second object in the shelf
contains a dictionary for the number of reads observed for each cell in the BAM file
that maps to the mitochondrial genome.

    Usage:
        python get_nucleotide_counts_for_mt_reads.py <input_bam> <output_file>

    Example:
python get_nucleotide_counts_for_mt_reads.py /media/kyle/Kyle_SSD_01/grimes_lab/data/Example_Files/example_mito_mapped_possorted_genome_male_first_1000.bam \
/media/kyle/Kyle_SSD_01/grimes_lab/analysis/2020/2020_04_17_tmporal_ordering_mitochondria_variant_alleles/output/mt_test_male_output_data

'''
import pysam
import sys
import shelve

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage: ")
        print("\tpython get_nucleotide_counts_for_mt_reads.py <input_bam> <output_file_path>")
        sys.exit()

    path_example = sys.argv[1]

    bam_file = pysam.AlignmentFile(path_example, "rb")

    base_counts = {}
    reads_per_cell = {}
    for read in bam_file.fetch():
        if read.has_tag('CB'):
            if read.reference_name == 'MT':
                if not read.get_tag('CB') in base_counts:
                    base_counts[read.get_tag('CB')] = {'A':{},'T':{},'C':{},'G':{}}
                    reads_per_cell[read.get_tag('CB')] = 1
                else:
                    reads_per_cell[read.get_tag('CB')] += 1

                base_position = read.reference_start
                for base in read.query_alignment_sequence:
                    if base in base_counts[read.get_tag('CB')]:
                        if base_position in base_counts[read.get_tag('CB')][base]:
                            base_counts[read.get_tag('CB')][base][base_position] += 1
                        else:
                            base_counts[read.get_tag('CB')][base][base_position] = 1

                    base_position += 1

    bam_file.close()

    shelf_file = shelve.open(sys.argv[2])
    shelf_file['base_counts'] = base_counts
    shelf_file['reads_per_cell'] = reads_per_cell
    shelf_file.close()

    print("Input bam file: \n\t%s\n\nWas successfully processed and saved to:\n\t%s" % (sys.argv[1], sys.argv[2]))
