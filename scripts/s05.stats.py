# generate stats for each processing steps

import pandas as pd

from collections import Counter
from functools import reduce

import gzip
import re
import os, sys


def count_reads(in_fq, n=4):

	# support either text or gzip file
	try:
		f = gzip.open(in_fq, "rt")
	except:
		f = open(in_fq, "rt")

	# count reads
	count = 0

	# fastq line
	fastq_n = 0

	for line in f:
		fastq_n += 1

		if fastq_n == n:
			count += 1
			fastq_n = 0

	# close file
	f.close()

	return count

def reads_length(in_fq, min_read_len, max_read_len, seq_n=2, n=4):
	# support either text or gzip file
	try:
		f = gzip.open(in_fq, "rt")
	except:
		f = open(in_fq, "rt")

	# read lengths
	reads_len = []

	# fastq line
	fastq_n = 0

	for line in f:
		fastq_n += 1

		if fastq_n == seq_n: # sequence line
			read_len = len(line.rstrip())

		if fastq_n == n:
			reads_len.append(read_len)
			fastq_n = 0

	# read length count
	reads_len_count = []

	counter = dict(Counter(reads_len))
	for i in range(min_read_len, max_read_len + 1):
		if i in counter:
			reads_len_count.append([i, counter[i]])
		else:
			reads_len_count.append([i, 0])

	f.close()

	return reads_len_count

def reads_length_null(min_read_len, max_read_len):
	reads_len_count = []

	for i in range(min_read_len, max_read_len + 1):
		reads_len_count.append([i, None])

	return reads_len_count

# simple parser for single-end star log file
def star_se_parser(star_log):
	input_reads = -1
	unique_reads = -1
	multi_reads = -1

	with open(star_log, "rt") as f:
		for line in f:
			if "Number of input reads" in line: # input reads
				m = re.search("^Number of input reads \|(.*)$", line.strip())
				if m:
					input_reads = int(m.group(1).strip())
			if "Uniquely mapped reads number" in line: # unique reads
				m = re.search("^Uniquely mapped reads number \|(.*)$", line.strip())
				if m:
					unique_reads = int(m.group(1).strip())
			if "Number of reads mapped to multiple loci" in line: # multi reads
				m = re.search("^Number of reads mapped to multiple loci \|(.*)$", line.strip())
				if m:
					multi_reads = int(m.group(1).strip())

	if any(count == -1 for count in [input_reads, unique_reads, multi_reads]):
		return [None, None, None]
	else:
		return [input_reads, unique_reads, multi_reads]


if __name__ == "__main__":
	# sample sheet
	samples = pd.read_table(snakemake.config["samples"]).set_index("sample_id", drop=False)

	# params
	min_read_len = snakemake.config["read_sel"]["min_len"]
	max_read_len = snakemake.config["read_sel"]["max_len"]

	# input files
	s00_input_fastq_files = snakemake.input.s00_input_fastq
	s01_trimming_fastq_files = snakemake.input.s01_trimming_fastq
	s02_read_sel_fastq_files = snakemake.input.s02_read_sel_fastq
	s03_ncrna_fastq_files = snakemake.input.s03_ncrna_fastq
	s03_mrna_fastq_files = snakemake.input.s03_mrna_fastq
	s04_star_log_files = snakemake.input.s04_star_log
	s04_mrna_transcriptome_fastq_files = snakemake.input.s04_mrna_transcriptome_fastq

	# output files
	s02_read_sel_stat_file = snakemake.output.s02_read_sel_stat
	s03_ncrna_align_stat_file = snakemake.output.s03_ncrna_align_stat
	s04_mrna_align_stat_file = snakemake.output.s04_mrna_align_stat
	s04_mrna_ncrna_read_len_stat_file = snakemake.output.s04_mrna_ncrna_read_len_stat

	# s02 read selection stat
	s02_read_sel_stat = []
	for i, sample_id in enumerate(samples["sample_id"].tolist()):
		s00_input_fastq_file = s00_input_fastq_files[i]
		s01_trimming_fastq_file = s01_trimming_fastq_files[i]
		s02_read_sel_fastq_file = s02_read_sel_fastq_files[i]

		# make sure file names in same order
		if all(sample_id in file_name for file_name in [s00_input_fastq_file, s01_trimming_fastq_file, s02_read_sel_fastq_file]):
			s02_read_sel_stat.append([sample_id, count_reads(s00_input_fastq_file), count_reads(s01_trimming_fastq_file), 
				count_reads(s02_read_sel_fastq_file)])
		else:
			sys.exit("File order inconsistent with sample sheet")
	s02_read_sel_stat = pd.DataFrame(s02_read_sel_stat, columns=["sample_id", "raw_reads", "trimming", "read_sel"])
	s02_read_sel_stat.to_csv(s02_read_sel_stat_file, sep="\t", header=True, index=False)

	# s03 ncrna alignment stat
	s03_ncrna_align_stat = []
	for i, sample_id in enumerate(samples["sample_id"].tolist()):
		s02_read_sel_fastq_file = s02_read_sel_fastq_files[i]
		s03_ncrna_fastq_file = s03_ncrna_fastq_files[i]
		s03_mrna_fastq_file = s03_mrna_fastq_files[i]

		if all(sample_id in file_name for file_name in [s02_read_sel_fastq_file, s03_ncrna_fastq_file, s03_mrna_fastq_file]):
			s03_ncrna_align_stat.append([sample_id, count_reads(s02_read_sel_fastq_file), count_reads(s03_ncrna_fastq_file), 
				count_reads(s03_mrna_fastq_file)])
		else:
			sys.exit("File order inconsistent with sample sheet")
	s03_ncrna_align_stat = pd.DataFrame(s03_ncrna_align_stat, columns=["sample_id", "read_sel", "ncrna", "non_ncrna"])
	s03_ncrna_align_stat["ncrna_pct"] = s03_ncrna_align_stat.apply(lambda row: row["ncrna"] / row["read_sel"], axis=1)
	s03_ncrna_align_stat["non_ncrna_pct"] = s03_ncrna_align_stat.apply(lambda row: row["non_ncrna"] / row["read_sel"], axis=1)
	s03_ncrna_align_stat.to_csv(s03_ncrna_align_stat_file, sep="\t", header=True, index=False)

	# s04 mrna alignment stat
	s04_mrna_align_stat = []
	for i, sample_id in enumerate(samples["sample_id"].tolist()):
		s04_star_log_file = s04_star_log_files[i]

		# make sure file names in same order
		if all(sample_id in file_name for file_name in [s04_star_log_file]):
			s04_mrna_align_stat.append([sample_id] + star_se_parser(s04_star_log_file))
		else:
			sys.exit("File order inconsistent with sample sheet")
	s04_mrna_align_stat = pd.DataFrame(s04_mrna_align_stat, columns=["sample_id", "non_ncrna", "unique", "multi"])
	s04_mrna_align_stat = s04_mrna_align_stat.fillna(value=pd.np.nan)
	s04_mrna_align_stat["unmapped"] = s04_mrna_align_stat.apply(lambda row: row["non_ncrna"] - row["unique"] - row["multi"], 
		axis=1)
	s04_mrna_align_stat["unique_pct"] = s04_mrna_align_stat.apply(lambda row: row["unique"] / row["non_ncrna"], axis=1)
	s04_mrna_align_stat["multi_pct"] = s04_mrna_align_stat.apply(lambda row: row["multi"] / row["non_ncrna"], axis=1)
	s04_mrna_align_stat["unmapped_pct"] = s04_mrna_align_stat.apply(lambda row: row["unmapped"] / row["non_ncrna"], axis=1)
	s04_mrna_align_stat.to_csv(s04_mrna_align_stat_file, sep="\t", header=True, index=False)

	# s04 mrna/ncrna read length stat
	s04_mrna_ncrna_read_len_stat = []
	for i, (sample_id, seq_type) in enumerate(zip(samples["sample_id"].tolist(), samples["seq_type"].tolist())):
		s03_ncrna_fastq_file = s03_ncrna_fastq_files[i]
		s04_mrna_transcriptome_fastq_file = s04_mrna_transcriptome_fastq_files[i]

		if all(sample_id in file_name for file_name in [s03_ncrna_fastq_file, s04_mrna_transcriptome_fastq_file]):
			if seq_type == "ribo":
				s03_ncrna_len_dist = reads_length(s03_ncrna_fastq_file, min_read_len, max_read_len)
				s04_mrna_len_dist = reads_length(s04_mrna_transcriptome_fastq_file, min_read_len, max_read_len)
				s03_ncrna_len_dist = pd.DataFrame(s03_ncrna_len_dist, columns=["read_len", "%s_ncrna" %(sample_id)])
				s04_mrna_len_dist = pd.DataFrame(s04_mrna_len_dist, columns=["read_len", "%s_non_ncrna" %(sample_id)])
				s04_mrna_ncrna_read_len_stat.append(pd.merge(s03_ncrna_len_dist, s04_mrna_len_dist, how="outer", on="read_len").
					fillna(value=pd.np.nan))
			elif seq_type == "rna":
				s03_ncrna_len_dist = reads_length_null(min_read_len, max_read_len)
				s04_mrna_len_dist = reads_length_null(min_read_len, max_read_len)
				s03_ncrna_len_dist = pd.DataFrame(s03_ncrna_len_dist, columns=["read_len", "%s_ncrna" %(sample_id)])
				s04_mrna_len_dist = pd.DataFrame(s04_mrna_len_dist, columns=["read_len", "%s_non_ncrna" %(sample_id)])
				s04_mrna_ncrna_read_len_stat.append(pd.merge(s03_ncrna_len_dist, s04_mrna_len_dist, how="outer", on="read_len").
					fillna(value=pd.np.nan))
			else:
				sys.exit("Invalid seq_type found: %s. Allowed values are [\"ribo\", \"rna\"]." %(seq_type))
		else:
			sys.exit("File order inconsistent with sample sheet")
	s04_mrna_ncrna_read_len_stat = reduce(lambda x, y: pd.merge(x, y, how="outer", on="read_len"), s04_mrna_ncrna_read_len_stat)
	s04_mrna_ncrna_read_len_stat.to_csv(s04_mrna_ncrna_read_len_stat_file, sep="\t", header=True, index=False, na_rep="NaN")
