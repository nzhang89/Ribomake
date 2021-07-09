import pandas as pd
from snakemake.utils import validate, min_version

import os, sys

# load config and sample sheets
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample_id", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")


# report
report: "reports/overview.rst"


# target rules
rule all:
	input:
		os.path.join(config["out_folder"], "data_processing.html")

# load rules
include: "rules/s01.trimming.smk"
include: "rules/s02.read_sel.smk"
include: "rules/s03.ncrna_align.smk"
include: "rules/s04.1.mrna_align.smk"
include: "rules/s04.2.mrna_fastq.smk"
include: "rules/s04.3.mrna_index.smk"
include: "rules/s05.stats.smk"
include: "rules/s06.report.smk"
