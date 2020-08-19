import shelve
import pandas as pd

path_shelf_file = '/data/salomonis2/LabFiles/Kyle/Analysis/2020_04_21_mito_af_lineage_tracing_human_hsc/output/'
path_output = '/data/salomonis2/LabFiles/Kyle/Analysis/2020_04_21_mito_af_lineage_tracing_human_hsc/output/human_hsc_h34_TA_'

shelf_file = shelve.open(path_shelf_file+'temp_mito_af_h34_TA')
reads_per_cell = pd.Series(shelf_file['reads_per_cell'])
base_counts = shelf_file['base_counts']
shelf_file.close()

# Filter out cells with fewer than 200 mapped reads
# coi = "cells of interest"
coi = pd.Series(reads_per_cell[reads_per_cell >= 200].index)

# Build 4 data frames:
#	af_a	: 	allele frequency of a
#	af_t	:	allele frequency of t
#	af_c	:	allele frequency of c
# 	af_g	:	allele frequency of g
# 	tot_a	:	total number of a
# 	tot_t	:	total number of t
# 	tot_c	:	total number of c
# 	tot_g	:	total number of g
af_a = {}
af_t = {}
af_c = {}
af_g = {}
tot_a = {}
tot_t = {}
tot_c = {}
tot_g = {}
tot_cov = {}
for cell in coi:
  tmp_df = pd.DataFrame(base_counts[cell])
  af_a[cell] = (tmp_df['A'] / tmp_df.sum(axis=1)).to_dict()
  af_t[cell] = (tmp_df['T'] / tmp_df.sum(axis=1)).to_dict()
  af_c[cell] = (tmp_df['C'] / tmp_df.sum(axis=1)).to_dict()
  af_g[cell] = (tmp_df['G'] / tmp_df.sum(axis=1)).to_dict()
  tot_a[cell] = base_counts[cell]['A']
  tot_t[cell] = base_counts[cell]['T']
  tot_c[cell] = base_counts[cell]['C']
  tot_g[cell] = base_counts[cell]['G']
  tot_cov[cell] = tmp_df.sum(axis=1).to_dict()
  
af_a_df = pd.DataFrame(af_a)
af_t_df = pd.DataFrame(af_t)
af_c_df = pd.DataFrame(af_c)
af_g_df = pd.DataFrame(af_g)
tot_a_df = pd.DataFrame(tot_a)
tot_t_df = pd.DataFrame(tot_t)
tot_c_df = pd.DataFrame(tot_c)
tot_g_df = pd.DataFrame(tot_g)

tot_cov_df = pd.DataFrame(tot_cov)
tot_cov_means = tot_cov_df.mean(axis=1)

# Select cells with at least 100 read mean coverage
selected_cells = list(tot_cov_means[tot_cov_means > 100].index)

tot_cov_means.to_csv(path_output+"mean_base_coverage.txt", sep="\t", index_label="UID")
af_a_df.iloc[selected_cells,:].to_csv(path_output+"allele_frequency_table_a.txt", sep="\t", index_label="UID")
af_t_df.iloc[selected_cells,:].to_csv(path_output+"allele_frequency_table_t.txt", sep="\t", index_label="UID")
af_c_df.iloc[selected_cells,:].to_csv(path_output+"allele_frequency_table_c.txt", sep="\t", index_label="UID")
af_g_df.iloc[selected_cells,:].to_csv(path_output+"allele_frequency_table_g.txt", sep="\t", index_label="UID")
tot_a_df.iloc[selected_cells,:].to_csv(path_output+"total_reads_table_a.txt", sep="\t", index_label="UID")
tot_t_df.iloc[selected_cells,:].to_csv(path_output+"total_reads_table_t.txt", sep="\t", index_label="UID")
tot_c_df.iloc[selected_cells,:].to_csv(path_output+"total_reads_table_c.txt", sep="\t", index_label="UID")
tot_g_df.iloc[selected_cells,:].to_csv(path_output+"total_reads_table_g.txt", sep="\t", index_label="UID")