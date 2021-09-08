import math
import pandas as pd


########################################################################################################################
##############################################   Parse outputs    ######################################################
########################################################################################################################

# Final outputs
outputs = expand("rlregions/rlregions_{group}{ext}", group=['All', 'S96', 'dRNH'], ext=[
  '.narrowPeak', '.bw', '.tsv', '_signal.tsv', ".bed"
])

########################################################################################################################
############################################   Helper Functions   ######################################################
########################################################################################################################

def get_minol(wildcards):
  """
  Finds the cutoff for peak calling based on number of input samples
  
  Chosen based on running macs3 bdgpeakcall with --cutoff-analysis
  """
  lns = pd.read_csv("rlregions/rlregion_manifest_" + wildcards.group + ".csv")
  res = math.ceil(lns.shape[0] * 0.20)
  if res < 5:
    res = 5
  return res


########################################################################################################################
##############################################   Main pipeline    ######################################################
#########################################################################################################################


rule output:
  input: outputs


rule clean_rlregions:
  input: 
    rl="rlregions/hg38_manifest_{group}.csv",
    np="rlregions/rlregions_{group}.narrowPeak"
  output: 
    tsv="rlregions/rlregions_{group}.tsv",
    bed="rlregions/rlregions_{group}.bed",
    sigtsv="rlregions/rlregions_{group}_signal.tsv"
  script: "scripts/finalizeRLRegions.R"


rule callpeak_from_bdg:
  input: "rlregions/counts.sorted_{group}.bdg"
  output: "rlregions/rlregions_{group}.narrowPeak"
  log: "logs/rlregions__callpeak_from_bdg_{group}.log"
  conda: "envs/macs.yaml"
  params:
    minol=get_minol
  shell: """
    macs3 bdgpeakcall -c {params.minol} -i {input} -o {output}
  """

rule toBigWig:
  input:
    bdg="rlregions/counts.sorted_{group}.bdg",
    chromsize="rlregions/hg38.chrom.sizes"
  output: "rlregions/rlregions_{group}.bw"
  log: "logs/rlregions__toBigWig_{group}.log"
  conda: "envs/bedgraph_to_bigwig.yaml"
  shell: """
    bedGraphToBigWig {input.bdg} {input.chromsize} {output}
  """

rule sortBed:
  input: "rlregions/counts.unsorted_{group}.bdg"
  output: "rlregions/counts.sorted_{group}.bdg"
  log: "logs/rlregions__sortBed_{group}.log"
  conda: "envs/bedgraph_to_bigwig.yaml"
  shell: """
    bedtools sort -i {input} > {output}
  """


rule toBedGraph:
  input: "rlregions/counts_{group}.bed"
  output: "rlregions/counts.unsorted_{group}.bdg"
  shell: """
    awk 'BEGIN {{FS="\t"; OFS="\t"}} {{if ($5 > 0) print $1, $2, $3, $5}}' {input} > {output}
  """


rule bedtools_intersect:
  output: "rlregions/counts_{group}.bed"
  input:
    manifest = "rlregions/rlregion_manifest_{group}.csv",
    windows = "rlregions/10bp_windows.bed"
  log: "logs/rlregions__bedtools_intersect_{group}.log"
  conda: "envs/bedtools.yaml"
  shell: """
    intersectBed -a {input.windows} -b $(tr '\n' ' ' < {input.manifest}) -c | \
    awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, $2, $3, $4, $5, "."}}' - > {output}
  """

rule makewindows:
  input: "rlregions/hg38.chrom.sizes"
  output: "rlregions/10bp_windows.bed"
  log: "logs/rlregions__makewindows_10bp.log"
  conda: "envs/bedtools.yaml"
  shell: "bedtools makewindows -g {input} -w 10 -i srcwinnum | sort -k 1,1 -k2,2n > {output}"


rule download_chromsizes:
  output: "rlregions/hg38.chrom.sizes"
  shell: "wget -O {output} ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
