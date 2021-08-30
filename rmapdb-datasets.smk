########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################


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
  "rmapdb-data/tables/rmap_rlfs.rda",
]

# Final outputs
outputs = UCSC_data + data_dumps + rmapdb_tables
            
########################################################################################################################
############################################   Helper Functions   ######################################################
########################################################################################################################



########################################################################################################################
##############################################   Main pipeline    ######################################################
#########################################################################################################################

rule output:
  input: outputs
  
  
rule rmap_rlfs:
  input: 
  output: "rmapdb-data/tables/rmap_rlfs.rda"





























