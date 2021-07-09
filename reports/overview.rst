This report summarizes the Ribosome profiling data processing Snakemake workflow, including adapter trimming, read length selection, ncRNA reads alignment, and non-ncRNA reads alignment.

All analysis results and log files are under the specified output folder: **{{snakemake.config["out_folder"]}}**, organized in the following structure:

* **s00.logs**: command execution log files for several steps;
* **s01.trimming**: adapter trimming results;
* **s02.read_sel**: read length selection results;
* **s03.ncrna_align**: ncRNA Bowtie2 alignment results;
* **s04.mrna_align**: non-ncRNA STAR alignment results;
* **s05.stats**: summary statistics for each processing step;
