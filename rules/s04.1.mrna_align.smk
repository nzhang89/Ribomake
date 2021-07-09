# mrna alignment step

rule s04_mrna_align_1:
	input:
		fastq=os.path.join(config["out_folder"], "s03.ncrna_align", "{sample}", "{sample}.trimmed.len_sel.mrna.fastq.gz")
	output:
		mrna_genome_bam=os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Aligned.sortedByCoord.out.bam"), 
		mrna_transcriptome_bam=os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Aligned.toTranscriptome.out.bam"),
		star_log=os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Log.final.out")
	params:
		main=config["mrna_align"]["params"], 
		index=config["mrna_align"]["index"], 
		gtf=config["mrna_align"]["gtf"], 
		out_prefix=os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", "{sample}.trimmed.len_sel.mrna.")
	threads:
		config["mrna_align"]["n_core"]
	conda:
		"../envs/all_env.yaml"
	log:
		os.path.join(config["out_folder"], "s00.logs", "s04.mrna_align", "{sample}.step1.log")
	shell:
		"STAR {params.main} --runThreadN {threads} --sjdbGTFfile {params.gtf} --genomeDir {params.index} "
		"--outFileNamePrefix {params.out_prefix} --readFilesIn {input} > {log} 2>&1"
