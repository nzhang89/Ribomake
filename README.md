# Ribomake: A Snakemake Workflow for Processing Ribosome Profiling Data

We present here Ribomake, a simple Snakemake workflow, for processing ribosome profiling data, including adapter trimming, non-coding RNA alignment, reference genome/transcriptome alignment, and data processing report. Ribomake is built and tested under 64-bit Linux system.

### Usage

#### Step 1: Install Snakemake

The first step is to install snakemake. We recommend using Conda to handle software dependencies. Please refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda) for more details on installing Snakemake via Conda.

Here, we show the basic installation process for snakemake 6.5.2 using Conda and Mamba.

```shell
# create conda environment named snakemake and install mamba
conda create -c conda-forge -n snakemake mamba

# activate snakemake environment
conda activate snakemake

# install snakemake into this environment
mamba install -c conda-forge -c bioconda snakemake=6.5.2

# make sure the following steps are also in this environment 
```

#### Step 2: Prepare input files

The Ribomake workflow requires some input files to be specified in either `SetupSnakemake.py` script or `config.yaml` file, listed below.

* Sample sheet file: a simple sample sheet file in TSV format (tab separated) containing four columns named as `sample_id`, `fq_path`, `seq_type`, and `adapter`. Columns `sample_id`, `fq_path`, and `seq_type` are required. Please **DO NOT** change the column names.

  * `sample_id`: stores unique ID/name for each sample;
  * `fq_path`: fastq file path for each sample;
  * `seq_type`: sequencing type. Possible values are "ribo" (for ribo-seq) or "rna" (for rna-seq);
  * `adapter`: specify here if some samples have different adapter sequences;

* Non-coding RNA Bowtie2 index: Ribo-seq reads might contain non-coding RNA contaminations (rRNA, tRNA, snRNA, snoRNA, etc.). We first align reads to these non-coding RNAs using Bowtie2. So you should prepare a Bowtie2 index file for the ncRNAs. ncRNA sequences for many species can be downloaded from Ensembl and tRNA sequences can be downloaded from specific tRNA databases.

* Reference genome STAR index: After aligning reads to ncRNAs, the remaining reads will be aligned to the reference genome using STAR. So you should also prepare a STAR index for the reference genome/transcriptome.

* Reference genome annotation GTF file: This should be the same GTF file used to generate STAR index.

#### Step 3: Set up Snakemake folder

First, clone this repository.

```shell
# clone and enter the repository
git clone https://github.com/nzhang89/Ribomake.git
cd Ribomake
```

Next, we will set up a Snakemake workflow folder. There are two ways to achieve this.

* Method 1 (Recommended):

Execute `SetupSnakemake.py` script. Run `python SetupSnakemake.py -h` to see the full list of arguments. The `snakemake_folder` is where the codes of the workflow are put into. It should be different from the current downloaded repository folder. The `analysis_folder` is where the output files and report should be put into.

```
usage: SetupSnakemake.py [-h] [--skip_adapter {False,True}] [--adapter ADAPTER] [--cutadapt_ncore CUTADAPT_NCORE]
                         [--cutadapt_params CUTADAPT_PARAMS] [--min_read_len MIN_READ_LEN] [--max_read_len MAX_READ_LEN]
                         --ncrna_index NCRNA_INDEX [--bowtie2_ncore BOWTIE2_NCORE] [--bowtie2_params BOWTIE2_PARAMS] --mrna_index
                         MRNA_INDEX --gtf GTF [--star_ncore STAR_NCORE] [--star_params STAR_PARAMS] [--python PYTHON]
                         [--pandas PANDAS] [--cutadapt CUTADAPT] [--bowtie2 BOWTIE2] [--star STAR] [--samtools SAMTOOLS]
                         [--shiny SHINY] [--rmarkdown RMARKDOWN] [--dt DT] [--plotly PLOTLY]
                         sample_sheet snakemake_folder analysis_folder

Set up Snakemake directory for ribosome profiling data analysis

positional arguments:
  sample_sheet                       Sample sheet file path (recommend using absolute path) (required)
  snakemake_folder                   Snakemake folder path (recommend using absolute path) (required)
  analysis_folder                    Data processing and analysis folder path (recommend using absolute path) (required)

optional arguments:
  -h, --help                         show this help message and exit

Adapter trimming step:
  --skip_adapter {False,True}        skip adapter trimming step (default: False)
  --adapter ADAPTER                  Ribosome profiling sequencing adapter for all samples (default: AGATCGGAAGAGCACACGTCT)
  --cutadapt_ncore CUTADAPT_NCORE    Number of CPU cores for each Cutadapt job (default: 4)
  --cutadapt_params CUTADAPT_PARAMS  Cutadapt trimming parameters (default: "-q 3 -e 0.1 --max-n 0.5 -m 15")

Read length selection step:
  --min_read_len MIN_READ_LEN        Minimum read length to keep (default: 25)
  --max_read_len MAX_READ_LEN        Maximum read length to keep (default: 35)

ncRNA alignment step:
  --ncrna_index NCRNA_INDEX          Path of Bowtie2 index for ncRNAs (required)
  --bowtie2_ncore BOWTIE2_NCORE      Number of CPU cores for each Bowtie2 job (default: 4)
  --bowtie2_params BOWTIE2_PARAMS    Bowtie2 parameters (default: "-k 1 --local --no-unal")

mRNA alignment step:
  --mrna_index MRNA_INDEX            Path of STAR index for transcriptome alignment (required)
  --gtf GTF                          Path of GTF annotation file used for STAR index (required)
  --star_ncore STAR_NCORE            Number of CPU cores for each STAR job (default: 4)
  --star_params STAR_PARAMS          STAR parameters. The default parameters are crucial for data processing and downstream analysis
                                     so please only add more options to the default if necessary (default: "--outSAMtype BAM
                                     SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --alignEndsType EndToEnd
                                     --readFilesCommand zcat")

Software package versions:
  --python PYTHON                    Python version (default: 3.6.8)
  --pandas PANDAS                    pandas version (default: 0.25.3)
  --cutadapt CUTADAPT                Cutadapt version (default: 2.6)
  --bowtie2 BOWTIE2                  Bowtie2 version (default: 2.4.4)
  --star STAR                        STAR version (default: 2.7.3a)
  --samtools SAMTOOLS                Samtools version (default: 1.12)
  --shiny SHINY                      shiny version (default: 1.4.0)
  --rmarkdown RMARKDOWN              rmarkdown version (default: 1.18)
  --dt DT                            dt version (default: 0.12)
  --plotly PLOTLY                    plotly version (default: 4.9.1)
```

* Method 2:
1. Copy the downloaded folder to your desired location.
2. Edit `config.yaml` and `envs/all_env.yaml` files. We provide templates for both files under `templates` folder.

After successfully setting up Snakemake folder (e.g. `/path/to/snakemake/folder`), switch to the Snakemake folder:
```shell
cd /path/to/snakemake/folder
```

Make sure the following files and folders exist and are organized in the correct structure.

```
.
├── config.yaml
├── envs
│   └── all_env.yaml
├── README.md
├── reports
│   ├── overview.rst
│   ├── s02_read_sel_stat.rst
│   ├── s03_ncrna_align_stat.rst
│   ├── s04_mrna_align_stat.rst
│   ├── s04_mrna_ncrna_read_len_stat.rst
│   └── s06_report.rst
├── rules
│   ├── s01.trimming.smk
│   ├── s02.read_sel.smk
│   ├── s03.ncrna_align.smk
│   ├── s04.1.mrna_align.smk
│   ├── s04.2.mrna_fastq.smk
│   ├── s04.3.mrna_index.smk
│   ├── s05.stats.smk
│   └── s06.report.smk
├── schemas
│   ├── config.schema.yaml
│   └── samples.schema.yaml
├── scripts
│   ├── s02.read_sel.py
│   ├── s05.stats.py
│   └── s06.report.Rmd
└── snakefile

5 directories, 23 files
```

#### Step 4: Run Snakemake workflow

We first create the conda environment with the specified software/packages and the corresponding versions. This step does not execute anything. It might take a few minutes to set up the conda environment. **Note** that for older versions of Snakemake, the command to create Conda environment might be different.

```shell
snakemake --use-conda --conda-create-envs-only --cores 1
```

Now we are ready to run Snakemake workflow locally:

```shell
snakemake --use-conda --cores 20
```

The `--cores` option controls the maximum number of CPU cores during workflow execution. For example, by default each STAR job uses 4 cores. If we set `--cores 20`, there are at most 5 STAR jobs running simultaneously.

If you plan to submit the workflow to clusters, please refer to [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

#### Step 5: Data processing and summary report

All analysis results and log files are under the specified analysis output folder (e.g. `/path/to/analysis/folder`), organized in the following structure. There is a data processing summary HTML report under the output folder.

| Folder              | Description                                   |
| ------------------- | --------------------------------------------- |
| **s00.logs**        | command execution log files for several steps |
| **s01.trimming**    | adapter trimming results                      |
| **s02.read_sel**    | read length selection results                 |
| **s03.ncrna_align** | ncRNA Bowtie2 alignment results               |
| **s04.mrna_align**  | non-ncRNA STAR alignment results              |
| **s05.stats**       | summary statistics for each processing step   |

Snakemake also allows us to generate a self-contained HTML report to our analysis output folder (e.g. `/path/to/analysis/folder`):
```shell
snakemake --report /path/to/analysis/folder/snakemake_report.html
```
The report contains the workflow summary, execution statistics, details for each step, and configuration information. In addition, we also embed the summary statistics for each step and each sample to the report under **Results** section. All files are downloadable tab-delimited text files.

#### Step 6: Routine data processing

If you plan to run the workflow on multiple datasets, you do not need to set up Snakemake folder every time. Instead, simply edit the `config.yaml` file (e.g. change sample sheet file, change output folder, etc.), save it in a different name, and run Snakemake using

```shell
snakemake --use-conda --cores 20 --configfile /path/to/new/config/file
```
