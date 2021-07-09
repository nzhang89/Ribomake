import gzip
import os, sys

# filter read length
# in_fq should be a text file (*.fastq or *.fq) or gzipped file (*.fastq.gz or *.fq.gz)
# out_fq is always a gzipped file (*.fastq.gz)
def filter_reads(in_fq, min_read_len, max_read_len, out_fq, chunk_size=10000000, n=4):
	# support either text or gzip file
	try:
		fi = gzip.open(in_fq, "rt")
	except:
		fi = open(in_fq, "rt")

	# output is always gzip
	fo = gzip.open(out_fq, "wt") # allow append

	# use chunk to minimize file write
	chunk = []

	# store one fastq read
	read = []

	for line in fi:
		read.append(line.rstrip())
		if len(read) == n:
			if len(read[1]) >= min_read_len and len(read[1]) <= max_read_len:
				chunk.append("\n".join(read))
			
			read.clear()

		if len(chunk) >= chunk_size:
			fo.write("\n".join(chunk))
			fo.write("\n")
			chunk.clear()
	fo.write("\n".join(chunk)) # write remaining reads

	# close files
	fi.close()
	fo.close()

def copy_reads(in_fq, out_fq, chunk_size=10000000, n=4):
	# support either text or gzip file
	try:
		fi = gzip.open(in_fq, "rt")
	except:
		fi = open(in_fq, "rt")

	# output is always gzip
	fo = gzip.open(out_fq, "wt") # allow append

	# use chunk to minimize file write
	chunk = []

	# store one fastq read
	read = []
	
	for line in fi:
		read.append(line.rstrip())
		if len(read) == n:
			chunk.append("\n".join(read))
			read.clear()

		if len(chunk) >= chunk_size:
			fo.write("\n".join(chunk))
			fo.write("\n")
			chunk.clear()
	fo.write("\n".join(chunk)) # write remaining reads

	# close files
	fi.close()
	fo.close()


if __name__ == "__main__":
	in_fq = snakemake.input.fastq
	out_fq = snakemake.output.fastq

	min_read_len = snakemake.params.min_len
	max_read_len = snakemake.params.max_len
	seq_type = snakemake.params.seq_type

	if seq_type == "ribo":
		filter_reads(in_fq, min_read_len, max_read_len, out_fq)
	elif seq_type == "rna":
		copy_reads(in_fq, out_fq)
	else:
		sys.exit("Invalid seq_type found: %s. Allowed values are [\"ribo\", \"rna\"]." %(seq_type))
