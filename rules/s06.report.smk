rule s06_report:
	input:
		s02_read_sel_stat=os.path.join(config["out_folder"], "s05.stats", "s02.read_sel.stats.txt"), 
		s03_ncrna_align_stat=os.path.join(config["out_folder"], "s05.stats", "s03.ncrna_align.stats.txt"), 
		s04_mrna_align_stat=os.path.join(config["out_folder"], "s05.stats", "s04.mrna_align.stats.txt"), 
		s04_mrna_ncrna_read_len_stat=os.path.join(config["out_folder"], "s05.stats", "s04.mrna_ncrna.read_len.stats.txt")
	output:
		s06_report=report(os.path.join(config["out_folder"], "data_processing.html"), 
			"../reports/s06_report.rst", category="5 Data Processing Report")
	conda:
		"../envs/all_env.yaml"
	script:
		"../scripts/s06.report.Rmd"
