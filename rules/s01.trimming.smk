# adapter trimming step

def get_fastq_s01(wildcards):
	return samples.loc[wildcards.sample, "fq_path"]

def get_adapter_s01(wildcards):
	adapter = config["trimming"]["adapter"]

	per_sample_adapter = samples.loc[wildcards.sample, "adapter"]
	if not pd.isnull(per_sample_adapter):
		adapter = per_sample_adapter

	return adapter


rule s01_trimming:
	input:
		fastq=get_fastq_s01
	output:
		fastq=os.path.join(config["out_folder"], "s01.trimming", "{sample}", "{sample}.trimmed.fastq.gz")
	params:
		main=config["trimming"]["params"], 
		adapter=lambda wildcards: get_adapter_s01(wildcards)
	threads:
		config["trimming"]["n_core"]
	conda:
		"../envs/all_env.yaml"
	log:
		os.path.join(config["out_folder"], "s00.logs", "s01.trimming", "{sample}.log")
	shell:
		"cutadapt {params.main} -j {threads} -a {params.adapter} -o {output.fastq} {input.fastq} > {log} 2>&1"
