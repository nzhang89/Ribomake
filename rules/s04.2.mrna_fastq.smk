# mrna alignment step

rule s04_mrna_align_2:
	input:
		mrna_transcriptome_bam=os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Aligned.toTranscriptome.out.bam")
	output:
		mrna_transcriptome_fastq=os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Aligned.toTranscriptome.out.fastq.gz")
	conda:
		"../envs/all_env.yaml"
	log:
		os.path.join(config["out_folder"], "s00.logs", "s04.mrna_align", "{sample}.step2.log")
	shell:
		"(samtools fastq {input.mrna_transcriptome_bam} | gzip -c > {output.mrna_transcriptome_fastq}) > {log} 2>&1"
