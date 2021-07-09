# ncrna alignment step

rule s03_ncrna_align:
	input:
		fastq=os.path.join(config["out_folder"], "s02.read_sel", "{sample}", "{sample}.trimmed.len_sel.fastq.gz")
	output:
		ncrna_fastq=os.path.join(config["out_folder"], "s03.ncrna_align", "{sample}", "{sample}.trimmed.len_sel.ncrna.fastq.gz"), 
		mrna_fastq=os.path.join(config["out_folder"], "s03.ncrna_align", "{sample}", "{sample}.trimmed.len_sel.mrna.fastq.gz"), 
		ncrna_bam=os.path.join(config["out_folder"], "s03.ncrna_align", "{sample}", "{sample}.trimmed.len_sel.ncrna.bam")
	params:
		main=config["ncrna_align"]["params"], 
		index=config["ncrna_align"]["index"], 
	threads:
		config["ncrna_align"]["n_core"]
	conda:
		"../envs/all_env.yaml"
	log:
		os.path.join(config["out_folder"], "s00.logs", "s03.ncrna_align", "{sample}.log")
	shell:
		"(bowtie2 {params.main} -p {threads} --al-gz {output.ncrna_fastq} --un-gz {output.mrna_fastq} -x {params.index} "
		"-U {input} | samtools view -b -o {output.ncrna_bam} -) > {log} 2>&1"
