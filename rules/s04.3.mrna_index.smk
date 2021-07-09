# mrna alignment step

rule s04_mrna_align_3:
	input:
		mrna_genome_bam = os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Aligned.sortedByCoord.out.bam")
	output:
		mrna_genome_bam_index=os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Aligned.sortedByCoord.out.bam.bai")
	conda:
		"../envs/all_env.yaml"
	log:
		os.path.join(config["out_folder"], "s00.logs", "s04.mrna_align", "{sample}.step3.log")
	shell:
		"(samtools index {input.mrna_genome_bam} {output.mrna_genome_bam_index}) > {log} 2>&1"
