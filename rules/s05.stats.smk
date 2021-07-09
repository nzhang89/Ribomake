# data processing report

def get_fastq_s05(skip_trimming):
	if skip_trimming:
		return samples["fq_path"]
	else:
		return expand(os.path.join(config["out_folder"], "s01.trimming", "{sample}", "{sample}.trimmed.fastq.gz"), 
			sample=samples["sample_id"])

rule s05_stats:
	input:
		s00_input_fastq=samples["fq_path"], 
		s01_trimming_fastq=get_fastq_s05(config["trimming"]["skip"]), 
		s02_read_sel_fastq=expand(os.path.join(config["out_folder"], "s02.read_sel", "{sample}", 
			"{sample}.trimmed.len_sel.fastq.gz"), sample=samples["sample_id"]), 
		s03_ncrna_fastq=expand(os.path.join(config["out_folder"], "s03.ncrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.ncrna.fastq.gz"), sample=samples["sample_id"]), 
		s03_mrna_fastq=expand(os.path.join(config["out_folder"], "s03.ncrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.fastq.gz"), sample=samples["sample_id"]), 
		s04_star_log=expand(os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Log.final.out"), sample=samples["sample_id"]), 
		s04_mrna_transcriptome_fastq=expand(os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Aligned.toTranscriptome.out.fastq.gz"), sample=samples["sample_id"]), 
		s04_mrna_genome_bam_index=expand(os.path.join(config["out_folder"], "s04.mrna_align", "{sample}", 
			"{sample}.trimmed.len_sel.mrna.Aligned.sortedByCoord.out.bam.bai"), sample=samples["sample_id"])
	output:
		s02_read_sel_stat=report(os.path.join(config["out_folder"], "s05.stats", "s02.read_sel.stats.txt"), 
			"../reports/s02_read_sel_stat.rst", category="1 Adapter Trimming and Read Length Selection Stats"), 
		s03_ncrna_align_stat=report(os.path.join(config["out_folder"], "s05.stats", "s03.ncrna_align.stats.txt"), 
			"../reports/s03_ncrna_align_stat.rst", category="2 ncRNA Alignment Stats"), 
		s04_mrna_align_stat=report(os.path.join(config["out_folder"], "s05.stats", "s04.mrna_align.stats.txt"), 
			"../reports/s04_mrna_align_stat.rst", category="3 mRNA alignment Stats"), 
		s04_mrna_ncrna_read_len_stat=report(os.path.join(config["out_folder"], "s05.stats", "s04.mrna_ncrna.read_len.stats.txt"), 
			"../reports/s04_mrna_ncrna_read_len_stat.rst", category="4 Read Length Distribution Stats")
	conda:
		"../envs/all_env.yaml"
	script:
		"../scripts/s05.stats.py"
