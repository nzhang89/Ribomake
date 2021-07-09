# read selection step

def get_fastq_s02(wildcards):
	if config["trimming"]["skip"]: # adapter trimming step is skipped
		return samples.loc[wildcards.sample, "fq_path"]
	else:
		return os.path.join(config["out_folder"], "s01.trimming", "%s" %(wildcards.sample), 
			"%s.trimmed.fastq.gz" %(wildcards.sample))

def get_seqtype_s02(wildcards):
	return samples.loc[wildcards.sample, "seq_type"]

rule s02_read_sel:
	input:
		fastq=get_fastq_s02
	output:
		fastq=os.path.join(config["out_folder"], "s02.read_sel", "{sample}", "{sample}.trimmed.len_sel.fastq.gz")
	params:
		min_len=config["read_sel"]["min_len"], 
		max_len=config["read_sel"]["max_len"], 
		seq_type=get_seqtype_s02
	conda:
		"../envs/all_env.yaml"
	script:
		"../scripts/s02.read_sel.py"
