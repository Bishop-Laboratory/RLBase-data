########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

import json
# config="rmap-data/rmap/config.json"
# config=json.load(open(config))

RSEQ_LOC="../RSeqCLI/"
RSEQR_LOC="../RSeqR/"
RTBL_LOC="rmapdb-data/tables/"
RSOUT_LOC= "rmapdb-data/rmap/"
DD_LOC="rmapdb-data/data-dumps/"
MISC_LOC="rmapdb-data/misc/"
CATALOG="RMapDB_manifest_27082021.xlsx"

sample=config['experiment']
genome=config['genome']

# Things that need to get made...
# From RSeq
"rmapdb-data/report/html/{sample}_{genome}.html"
"rmapdb-data/report/data/{sample}_{genome}.rda"
"rmapdb-data/peaks/{sample}_{genome}.broadPeak"
"rmapdb-data/coverage/{sample}_{genome}.bw"
"rmapdb-data/bam/{sample}_{genome}.bam"
"rmapdb-data/bam/{sample}_{genome}.bam.bai"
# Expression
"rmapdb-data/expression/{sample}_{genome}/quant.sf"
# RLFS
"rmapdb-data/rlfs-beds/{genome}.rlfs.bed"


########################################################################################################################
##############################################   Parse outputs    ######################################################
########################################################################################################################

# UCSC Genome Browser
UCSC_data = [
  "rloop_regions.bw",
  "rloop_regions.bed"
]
UCSC_data = [MISC_LOC + ud for ud in UCSC_data]


# Data dumps
data_dumps = [
  "rmapdb_peaks.tar.xz",
  "rmapdb_coverage.tar.xz",
  "rmapdb_report_html.tar.xz",
  "rmapdb_report_data.tar.xz",
  "rmapdb_expression.tar.xz"
]
data_dumps = [DD_LOC + dd for dd in data_dumps]

# RMapDB Tables
rmapdb_tables = [
  # "rloop_regions.tsv.xz",
  # "rloop_signal.tsv.xz",
  # "gene_rl_overlap.tsv.xz",
  # "gene_expression.tsv.xz",
  # "gene_exp_to_rmap_samples.tsv.xz",
  "gene_exp_samples.tsv.xz",
  "rmap_samps.tsv.xz",
  "rmap_corr.tsv.xz",
  "rmap_pca.tsv.xz",
  "rmap_anno.tsv.xz",
  "rmap_rlfs.tar.xz"
]
rmapdb_tables = [RTBL_LOC + tbl for tbl in rmapdb_tables]

# Final outputs
outputs = UCSC_data + data_dumps + rmapdb_tables
outputs=rmapdb_tables
            
########################################################################################################################
############################################   Helper Functions   ######################################################
########################################################################################################################



########################################################################################################################
##############################################   Main pipeline    ######################################################
#########################################################################################################################

# Define sample-level inputs
rlfs_rdas=[RTBL_LOC + "rmap_rlfs/" + samp + "_" + gen + ".rda" for samp, gen in zip(sample, genome)]
peaks=[RSOUT_LOC + "peaks/" + samp + "_" + gen + ".broadPeak" for samp, gen in zip(sample, genome)]
coverage=[RSOUT_LOC + "coverage/" + samp + "_" + gen + ".bw" for samp, gen in zip(sample, genome)]


rlfs_rdas=rlfs_rdas[1:5]
peaks=peaks[1:5]
coverage=coverage[1:5]

rule output:
  input: outputs
  
rule tarball:
  """Will tar any directory with the 'tarball.txt' file"""
  input: "{path}/tarball.txt"
  output: "{path}.tar.xz"
  params:
    dir="{path}"
  shell: "tar cfJ {output} {params.dir}"

rule xz:
  """Will xz any file"""
  input: "{path}"
  output: "{path}.xz"
  params:
    dir="{path}"
  shell: "xz -f {input}"
  
rule get_rmap_samps:
  input: 
    manifest="rmap-data/rmap_manifest.csv",
    model=MISC_LOC + "fft_model.rda"
  output: RTBL_LOC + "rmap_samps.tsv"
  script: "scripts/compileRMapSamps.R"
  
rule build_model:
  input:
    manifest="rmap-data/rmap_manifest.csv",
    blacklist=MISC_LOC + "rmap_blacklist.csv",
    rlfs=rlfs_rdas
  output: MISC_LOC + "fft_model.rda"
  script: "scripts/doBuildModel.R"
  

rule make_blacklist:
  input:
    rmap="rmap-data/rmap_manifest.csv",
    rlfs=rlfs_rdas
  output: MISC_LOC + "rmap_blacklist.csv"
  conda: "envs/jupyteR.yaml"
  notebook: "notebooks/blacklist.r.ipynb"
  

rule get_pca:
  input: RTBL_LOC + "rmap_corr.tsv"
  output: RTBL_LOC + "rmap_pca.tsv"
  script: "scripts/doPCAAnalysis.R"

rule get_corr:
  input: peaks
  output: RTBL_LOC + "rmap_corr.tsv"
  script: "scripts/doCorrAnalysis.R"

rule gather_anno:
  input:  peaks
  output: RTBL_LOC + "rmap_anno.tsv"
  script: "scripts/doAnnoAnalysis.R"

rule prep_rlfs_out:
  """Aggregator for RLFS ahead of tarball"""
  input: rlfs_rdas
  output: RTBL_LOC + "rmap_rlfs/tarball.txt"
  shell: "touch {output}"

rule gather_rlfs:
  input: RSOUT_LOC+"peaks/{sample}_{genome}.broadPeak"
  output: RTBL_LOC+"rmap_rlfs/{sample}_{genome}.rda"
  script: "scripts/doAnalyzeRLFS.R"

rule rseq_run_norep:
  """RSeqCLI run with --no-report flag"""
  input: RSOUT_LOC+"config.json"
  output: 
    peak=RSOUT_LOC+"peaks/{sample}_{genome}.broadPeak",
    coverage=RSOUT_LOC+"coverage/{sample}_{genome}.bw",
    bam=RSOUT_LOC+"bam/{sample}_{genome}.bam",
    bai=RSOUT_LOC+"bam/{sample}_{genome}.bam.bai"
  shell: "RSeqCLI run rmap-data/rmap/ --bwamem2 -t 44 --no-report"
  
rule rseq_check_norep:
  """RSeqCLI run with --no-report flag"""
  input: RSOUT_LOC+"config.json"
  output: RSOUT_LOC+"checknorep.txt"
  shell: """
    RSeqCLI check rmap-data/rmap/ --bwamem2 -t 44 --no-report
    touch {output}
  """

rule rseq_build:
  """RSeqCLI build"""
  input: "rmap-data/rmap_manifest.csv"
  output: RSOUT_LOC+"config.json"
  shell: """
    RSeqCLI build rmap-data/rmap/ {input}
  """

rule makeManifests:
  input: CATALOG
  output: 
    rmap="rmap-data/rmap_manifest.csv"
  script: "scripts/doMakeManifest.R"


