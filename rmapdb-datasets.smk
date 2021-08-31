########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

RSEQ_LOC="../RSeqCLI/"
RSEQR_LOC="../RSeqR/"
CATALOG="RMapDB_manifest_27082021.xlsx"

sample=config['experiment']
genome=config['genome']
prep_rlfs_out_inputs=["rmapdb-data/tables/rmap_rlfs/" + samp + "_" + gen + ".rda" for samp,gen in zip(sample, genome)]
prep_rlfs_out_inputs=prep_rlfs_out_inputs[1:10]


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
  "rmapdb-data/misc/rloop_regions.bw",
  "rmapdb-data/misc/rloop_regions.bed"
]

# Data dumps
data_dumps = [
  "rmapdb-data/data-dumps/rmapdb_peaks.tar.xz",
  "rmapdb-data/data-dumps/rmapdb_coverage.tar.xz",
  "rmapdb-data/data-dumps/rmapdb_report_html.tar.xz",
  "rmapdb-data/data-dumps/rmapdb_report_data.tar.xz",
  "rmapdb-data/data-dumps/rmapdb_expression.tar.xz"
]

# RMapDB Tables
rmapdb_tables = [
  "rmapdb-data/tables/rloop_regions.tsv.xz",
  "rmapdb-data/tables/rloop_signal.tsv.xz",
  "rmapdb-data/tables/gene_rl_overlap.tsv.xz",
  "rmapdb-data/tables/gene_expression.tsv.xz",
  "rmapdb-data/tables/gene_exp_to_rmap_samples.tsv.xz",
  "rmapdb-data/tables/gene_exp_samples.tsv.xz",
  "rmapdb-data/tables/rmap_samps.tsv.xz",
  "rmapdb-data/tables/rmap_corr.tsv.xz",
  "rmapdb-data/tables/rmap_pca.tsv.xz",
  "rmapdb-data/tables/rmap_anno.tsv.xz",
  "rmapdb-data/tables/rmap_rlfs.tar.xz",
]

# Final outputs
outputs = UCSC_data + data_dumps + rmapdb_tables
            
########################################################################################################################
############################################   Helper Functions   ######################################################
########################################################################################################################



########################################################################################################################
##############################################   Main pipeline    ######################################################
#########################################################################################################################

outputs = "rmapdb-data/tables/rmap_rlfs.tar.xz"

rule output:
  input: outputs
  
rule tarball:
  """Will tar any directory with the 'tarball.txt' file"""
  input: "rmapdb-data/{path}/tarball.txt"
  output: "rmapdb-data/{path}.tar.xz"
  params:
    dir="rmapdb-data/{path}"
  shell: "tar cfJ {output} {params.dir}"
  
rule prep_rlfs_out:
  """Aggregator for RLFS ahead of tarball"""
  input: prep_rlfs_out_inputs
  output: "rmapdb-data/tables/rmap_rlfs/tarball.txt"
  shell: "touch {output}"

rule gather_rlfs:
  input: "rmapdb-data/rmap/peaks/{sample}_{genome}.broadPeak"
  output: "rmapdb-data/tables/rmap_rlfs/{sample}_{genome}.rda"
  script: "scripts/doAnalyzeRLFS.R"

rule rseq_run_norep:
  """RSeqCLI run with --no-report flag"""
  input: "rmap-data/rmap/config.json"
  output: 
    peak="rmapdb-data/rmap/peaks/{sample}_{genome}.broadPeak",
    coverage="rmapdb-data/rmap/coverage/{sample}_{genome}.bw",
    bam="rmapdb-data/rmap/bam/{sample}_{genome}.bam",
    bai="rmapdb-data/rmap/bam/{sample}_{genome}.bam.bai"
  shell: "RSeqCLI run rmap-data/rmap/ --bwamem2 -t 44 --no-report"
  
rule rseq_check_norep:
  """RSeqCLI run with --no-report flag"""
  input: "rmap-data/rmap/config.json"
  output: "rmap-data/rmap/checknorep.txt"
  shell: """
    RSeqCLI check rmap-data/rmap/ --bwamem2 -t 44 --no-report
    touch {output}
  """

rule rseq_build:
  """RSeqCLI build"""
  input: "rmap-data/rmap_manifest.csv"
  output: "rmap-data/rmap/config.json"
  shell: """
    RSeqCLI build rmap-data/rmap/ {input}
  """

rule makeManifests:
  input: CATALOG
  output: 
    rmap="rmap-data/rmap_manifest.csv"
  script: "scripts/doMakeManifest.R"


