# set up snakemake directory

import argparse

from distutils.util import strtobool
from distutils.dir_util import copy_tree, remove_tree
from distutils.file_util import copy_file

import yaml

import os, sys, time


def load_default_params():
	return {
		"skip_adapter" : "False", 
		"adapter" : "AGATCGGAAGAGCACACGTCT", 
		"cutadapt_ncore": 4, 
		"cutadapt_params" : "-q 3 -e 0.1 --max-n 0.5 -m 15", 

		"min_read_len" : 25, 
		"max_read_len" : 35, 

		"bowtie2_ncore" : 4, 
		"bowtie2_params" : "-k 1 --local --no-unal", 

		"star_ncore" : 4, 
		"star_params" : " ".join(["--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts", 
			"--alignEndsType EndToEnd --readFilesCommand zcat"]), 

		"python" : "3.6.8", 
		"pandas" : "0.25.3",
		"cutadapt" : "2.6", 
		"bowtie2" : "2.4.4", 
		"star" : "2.7.3a", 
		"samtools" : "1.12", 
		"shiny": "1.4.0", 
		"rmarkdown" : "1.18", 
		"dt": "0.12", 
		"plotly": "4.9.1"
	}


def now():
	return time.strftime("%Y-%m-%d %H:%M:%S %Z", time.localtime())

def check_files(curr_dir):
	complete_rulefiles = ["s01.trimming.smk", "s02.read_sel.smk", "s03.ncrna_align.smk", "s04.1.mrna_align.smk", 
		"s04.2.mrna_fastq.smk", "s05.stats.smk"]
	complete_schemafiles = ["config.schema.yaml", "samples.schema.yaml"]
	complete_scriptfiles = ["s02.read_sel.py", "s05.stats.py"]

	configfile_check = os.path.exists(os.path.join(curr_dir, "config.yaml"))
	snakefile_check = os.path.exists(os.path.join(curr_dir, "snakefile"))
	envfile_check = os.path.exists(os.path.join(curr_dir, "envs", "all_env.yaml"))
	rulefiles_check = all(os.path.exists(os.path.join(curr_dir, "rules", rulefile)) for rulefile in complete_rulefiles)
	schemafiles_check = all(os.path.exists(os.path.join(curr_dir, "schemas", schemafile)) for schemafile in complete_schemafiles)
	scriptfiles_check = all(os.path.exists(os.path.join(curr_dir, "scripts", scriptfile)) for scriptfile in complete_scriptfiles)

	if not all([configfile_check, snakefile_check, envfile_check, rulefiles_check, schemafiles_check, scriptfiles_check]):
		print("", file=sys.stderr)

	if not configfile_check:
		sys.exit("File not exists: config.yaml under %s" %(curr_dir))
	if not snakefile_check:
		sys.exit("File not exists: snakefile under %s" %(curr_dir))
	if not envfile_check:
		sys.exit("File not exists: all_env.yaml under %s" %(os.path.join(curr_dir, "envs")))
	if not rulefiles_check:
		sys.exit("Snakemake rule files incomplete under %s. Complete rule files are [%s]" 
			%(os.path.join(curr_dir, "rules"), ", ".join(complete_rulefiles)))
	if not schemafiles_check:
		sys.exit("Snakemake schema files incomplete under %s. Complete schema files are [%s]" 
			%(os.path.join(curr_dir, "schemas"), ", ".join(complete_schemafiles)))
	if not scriptfiles_check:
		sys.exit("Snakemake script files incomplete under %s. Complete script files are [%s]" 
			%(os.path.join(curr_dir, "scripts"), ", ".join(complete_scriptfiles)))

def process_args():
	# default params
	default = load_default_params()

	# format arguments
	formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=132)
	parser = argparse.ArgumentParser(description="Set up Snakemake directory for ribosome profiling data analysis", 
		formatter_class=formatter)

	# arguments
	parser.add_argument("sample_sheet", type=str, 
		help="Sample sheet file path (recommend using absolute path) (required)")
	parser.add_argument("snakemake_folder", type=str, 
		help="Snakemake folder path (recommend using absolute path) (required)")
	parser.add_argument("analysis_folder", type=str, 
		help="Data processing and analysis folder path (recommend using absolute path) (required)")

	# trimming step
	s01_trimming = parser.add_argument_group("Adapter trimming step")
	s01_trimming.add_argument("--skip_adapter", type=str, default=default["skip_adapter"], choices=["False", "True"], 
		help="skip adapter trimming step (default: %s)" %(default["skip_adapter"]))
	s01_trimming.add_argument("--adapter", type=str, default=default["adapter"], 
		help="Ribosome profiling sequencing adapter for all samples (default: %s)" %(default["adapter"]))
	s01_trimming.add_argument("--cutadapt_ncore", type=int, default=default["cutadapt_ncore"], 
		help="Number of CPU cores for each Cutadapt job (default: %d)" %(default["cutadapt_ncore"]))
	s01_trimming.add_argument("--cutadapt_params", type=str, default=default["cutadapt_params"], 
		help="Cutadapt trimming parameters (default: \"%s\")" %(default["cutadapt_params"]))

	# read length selection step
	s02_read_sel = parser.add_argument_group("Read length selection step")
	s02_read_sel.add_argument("--min_read_len", type=int, default=default["min_read_len"], 
		help="Minimum read length to keep (default: %d)" %(default["min_read_len"]))
	s02_read_sel.add_argument("--max_read_len", type=int, default=default["max_read_len"], 
		help="Maximum read length to keep (default: %d)" %(default["max_read_len"]))
	
	# ncrna alignment step
	s03_ncrna_align = parser.add_argument_group("ncRNA alignment step")
	s03_ncrna_align.add_argument("--ncrna_index", type=str, required=True, 
		help="Path of Bowtie2 index for ncRNAs (required)")
	s03_ncrna_align.add_argument("--bowtie2_ncore", type=int, default=default["bowtie2_ncore"], 
		help="Number of CPU cores for each Bowtie2 job (default: %d)" %(default["bowtie2_ncore"]))
	s03_ncrna_align.add_argument("--bowtie2_params", type=str, default=default["bowtie2_params"], 
		help="Bowtie2 parameters (default: \"%s\")" %(default["bowtie2_params"]))

	# mrna alignment step
	s04_mrna_align = parser.add_argument_group("mRNA alignment step")
	s04_mrna_align.add_argument("--mrna_index", type=str, required=True, 
		help="Path of STAR index for transcriptome alignment (required)")
	s04_mrna_align.add_argument("--gtf", type=str, required=True, 
		help="Path of GTF annotation file used for STAR index (required)")
	s04_mrna_align.add_argument("--star_ncore", type=int, default=default["star_ncore"], 
		help="Number of CPU cores for each STAR job (default: %d)" %(default["star_ncore"]))
	s04_mrna_align.add_argument("--star_params", type=str, default=default["star_params"], 
		help=" ".join(["STAR parameters. The default parameters are crucial for data processing and downstream analysis", 
		"so please only add more options to the default if necessary (default: \"%s\")" %(default["star_params"])]))

	# software package version
	s00_package = parser.add_argument_group("Software package versions")
	s00_package.add_argument("--python", type=str, default=default["python"], 
		help="Python version (default: %s)" %(default["python"]))
	s00_package.add_argument("--pandas", type=str, default=default["pandas"], 
		help="pandas version (default: %s)" %(default["pandas"]))
	s00_package.add_argument("--cutadapt", type=str, default=default["cutadapt"], 
		help="Cutadapt version (default: %s)" %(default["cutadapt"]))
	s00_package.add_argument("--bowtie2", type=str, default=default["bowtie2"], 
		help="Bowtie2 version (default: %s)" %(default["bowtie2"]))
	s00_package.add_argument("--star", type=str, default=default["star"], 
		help="STAR version (default: %s)" %(default["star"]))
	s00_package.add_argument("--samtools", type=str, default=default["samtools"], 
		help="Samtools version (default: %s)" %(default["samtools"]))
	s00_package.add_argument("--shiny", type=str, default=default["shiny"], 
		help="shiny version (default: %s)" %(default["shiny"])), 
	s00_package.add_argument("--rmarkdown", type=str, default=default["rmarkdown"], 
		help="rmarkdown version (default: %s)" %(default["rmarkdown"])), 
	s00_package.add_argument("--dt", type=str, default=default["dt"], 
		help="dt version (default: %s)" %(default["dt"])), 
	s00_package.add_argument("--plotly", type=str, default=default["plotly"], 
		help="plotly version (default: %s)" %(default["plotly"]))

	args = parser.parse_args()

	return args

def setup_config(args, config_file, out_file):
	with open(config_file, "rt") as f:
		config = yaml.safe_load(f)

		config["out_folder"] = os.path.abspath(args.analysis_folder)
		config["samples"] = os.path.abspath(args.sample_sheet)

		config["trimming"]["skip"] = bool(strtobool(args.skip_adapter))
		config["trimming"]["adapter"] = args.adapter
		config["trimming"]["ncore"] = args.cutadapt_ncore
		config["trimming"]["params"] = args.cutadapt_params

		config["read_sel"]["min_len"] = args.min_read_len
		config["read_sel"]["max_len"] = args.max_read_len

		config["ncrna_align"]["index"] = os.path.abspath(args.ncrna_index)
		config["ncrna_align"]["n_core"] = args.bowtie2_ncore
		config["ncrna_align"]["params"] = args.bowtie2_params

		config["mrna_align"]["index"] = os.path.abspath(args.mrna_index)
		config["mrna_align"]["gtf"] = os.path.abspath(args.gtf)
		config["mrna_align"]["n_core"] = args.star_ncore
		config["mrna_align"]["params"] = args.star_params

		config["package"]["python"] = args.python
		config["package"]["pandas"] = args.pandas
		config["package"]["cutadapt"] = args.cutadapt
		config["package"]["bowtie2"] = args.bowtie2
		config["package"]["star"] = args.star
		config["package"]["samtools"] = args.samtools
		config["package"]["shiny"] = args.shiny
		config["package"]["rmarkdown"] = args.rmarkdown
		config["package"]["dt"] = args.dt
		config["package"]["plotly"] = args.plotly

	with open(out_file, "wt") as f:
		yaml.dump(config, f, default_flow_style=False, sort_keys=False, indent=2)

def setup_env(args, env_file, out_file):
	with open(env_file, "rt") as f:
		env = yaml.safe_load(f)

		package_ver = []
		package_ver.append("python=%s" %(args.python))
		package_ver.append("pandas=%s" %(args.pandas))
		package_ver.append("cutadapt=%s" %(args.cutadapt))
		package_ver.append("bowtie2=%s" %(args.bowtie2))
		package_ver.append("star=%s" %(args.star))
		package_ver.append("samtools=%s" %(args.samtools))
		package_ver.append("r-shiny=%s" %(args.shiny))
		package_ver.append("r-rmarkdown=%s" %(args.rmarkdown))
		package_ver.append("r-dt=%s" %(args.dt))
		package_ver.append("r-plotly=%s" %(args.plotly))

		env["dependencies"] = package_ver

	with open(out_file, "wt") as f:
		yaml.safe_dump(env, f, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
	# current folder
	curr_dir = os.path.dirname(os.path.realpath(__file__))
	
	# args
	args = process_args()
	check_files(curr_dir)

	# folder setup
	print("[%s] Setting up Snakemake and data analysis folder" %(now()), file=sys.stderr)

	# snakemake folder
	snakemake_folder = args.snakemake_folder
	if not os.path.exists(snakemake_folder): # create snakemake folder if not exists
		os.makedirs(snakemake_folder)
	else:
		print("[%s] Snakemake folder %s already exists. Snakemake related files will be overwritten." 
			%(now(), snakemake_folder), file=sys.stderr)
		if os.path.exists(os.path.join(snakemake_folder, "envs")):
			remove_tree(os.path.join(snakemake_folder, "envs"))
		if os.path.exists(os.path.join(snakemake_folder, "reports")):
			remove_tree(os.path.join(snakemake_folder, "reports"))
		if os.path.exists(os.path.join(snakemake_folder, "rules")):
			remove_tree(os.path.join(snakemake_folder, "rules"))
		if os.path.exists(os.path.join(snakemake_folder, "schemas")):
			remove_tree(os.path.join(snakemake_folder, "schemas"))
		if os.path.exists(os.path.join(snakemake_folder, "scripts")):
			remove_tree(os.path.join(snakemake_folder, "scripts"))
		if os.path.exists(os.path.join(snakemake_folder, "snakefile")):
			os.remove(os.path.join(snakemake_folder, "snakefile"))
		if os.path.exists(os.path.join(snakemake_folder, "config.yaml")):
			os.remove(os.path.join(snakemake_folder, "config.yaml"))
		if os.path.exists(os.path.join(snakemake_folder, ".snakemake")):
			remove_tree(os.path.join(snakemake_folder, ".snakemake"))

	# copy folders/files to snakemake folder
	if not os.path.exists(os.path.join(snakemake_folder, "envs")):
		os.makedirs(os.path.join(snakemake_folder, "envs"))
	setup_env(args, os.path.join(curr_dir, "envs", "all_env.yaml"), os.path.join(snakemake_folder, "envs", "all_env.yaml"))

	copy_tree(os.path.join(curr_dir, "reports"), os.path.join(snakemake_folder, "reports"))
	copy_tree(os.path.join(curr_dir, "rules"), os.path.join(snakemake_folder, "rules"))
	copy_tree(os.path.join(curr_dir, "schemas"), os.path.join(snakemake_folder, "schemas"))
	copy_tree(os.path.join(curr_dir, "scripts"), os.path.join(snakemake_folder, "scripts"))
	copy_file(os.path.join(curr_dir, "README.md"), snakemake_folder)
	copy_file(os.path.join(curr_dir, "snakefile"), snakemake_folder)

	setup_config(args, os.path.join(curr_dir, "config.yaml"), os.path.join(snakemake_folder, "config.yaml"))
	
	# data analysis folder
	analysis_folder = args.analysis_folder
	if not os.path.exists(analysis_folder):
		os.makedirs(analysis_folder)

	print("[%s] Finished Snakemake/data analysis folder setup. Switch to Snakemake folder and refer to README for running jobs" 
		%(now()), file=sys.stderr)
