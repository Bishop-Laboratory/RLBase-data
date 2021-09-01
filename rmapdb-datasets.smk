########################################################################################################################
##############################################   Parse inputs    #######################################################
########################################################################################################################

import json

# config="rmap-data/rmap/config.json"
# config=json.load(open(config))
RSEQ_LOC="../RSeqCLI/"
RSEQR_LOC="../RSeqR/"
CATALOG="RMapDB_manifest_27082021.xlsx"
CATALOG="RMapDB_manifest_27082021.small.xlsx"

sample=config['experiment']
genome=config['genome']
mode=config['mode']

########################################################################################################################
##############################################   Parse outputs    ######################################################
########################################################################################################################

# UCSC Genome Browser
other_data = [
  "rlregions.bw",
  "rlregions.bed",
  "rmap_rlfs.tar.xz"
]
other_data = ["misc/" + ud for ud in other_data]


# Data dumps
data_dumps = [
  "rmapdb_peaks",
  "rmapdb_coverage",
  "rmapdb_report_html",
  "rmapdb_report_data",
  "rmapdb_expression"
]
data_dumps = ["data-dumps/" + dd + ".tar.xz" for dd in data_dumps]

# RMapDB Tables
rmapdb_tables = [
  "rlregions_signal.tsv",
  "rlregions.tsv",
  "gene_expression.tsv",
  "rmap_samps.tsv",
  "rmap_corr.tsv",
  "rmap_pca.tsv",
  "rmap_anno.tsv"
]
rmapdb_tables = ["tables/" + tbl for tbl in rmapdb_tables]

# Final outputs
outputs = other_data + rmapdb_tables #+ data_dumps
            
########################################################################################################################
############################################   Helper Functions   ######################################################
########################################################################################################################



########################################################################################################################
##############################################   Main pipeline    ######################################################
#########################################################################################################################

# Define sample-level inputs
rlfs_rdas=["tables/rmap_rlfs/" + elem + "_" + genome[idx] + ".rda" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
peaks=["rmap/peaks/" + elem + "_" + genome[idx] + ".broadPeak" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
bams=["rmap/bam/" + elem + "_" + genome[idx] + "/" + elem + "_" + genome[idx] + ".bam" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
coverage=["rmap/coverage/" + elem + "_" + genome[idx] + ".bw" for idx, elem in enumerate(sample) if mode[idx] not in ["RNA-Seq", "RNA-seq"]]
quant=["rmap/quant/" + 'quant/' + elem + "_" + genome[idx] + "/quant.sf" for idx, elem in enumerate(sample) if mode[idx] in ["RNA-Seq", "RNA-seq"]]


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
  
  
rule clean_rlregion_signal:
  input: "rlregions/rlregions__signal_raw.tsv"
  output: "tables/rlregions_signal.tsv"
  script: "scripts/cleanRLRegionsSignal.R"

rule clean_rlregions:
  input: 
    rl="rlregions/rlregions__signal_raw.tsv",
    np="rlregions/rlregions.narrowPeak"
  output: 
    tsv="tables/rlregions.tsv",
    bed="misc/rlregions.bed"

  
rule rlregion_signal:
  input:
    signals=coverage,
    peaks="rlregions/rlregions.narrowPeak"
  output:
    npz="rlregions/rlregions_signal.npz",
    tsv="rlregions/rlregions__signal_raw.tsv"
  conda: "envs/deeptools.yaml"
  threads: 30
  log: "logs/rlregions__deeptools.log"
  shell: """
  (
    multiBigwigSummary BED-file -b {input.signals} --BED {input.peaks} \
    -p {threads} -v -o {output.npz} --outRawCounts {output.tsv}
  ) &> {log}
  """
  
  
rule callpeak_from_bdg:
  input: "rlregions/counts.sorted.bdg"
  output: "rlregions/rlregions.narrowPeak"
  log: "logs/rlregions__callpeak_from_bdg.log"
  conda: "envs/macs.yaml"
  shell: """
    macs3 bdgpeakcall -i {input} -o {output}
  """
  

rule toBigWig:
  input:
    bdg="rlregions/counts.sorted.bdg",
    chromsize="rlregions/hg38.chrom.sizes"
  output: "misc/rlregions.bw"
  log: "logs/rlregions__toBigWig.log"
  conda: "envs/bedgraph_to_bigwig.yaml"
  shell: """
    bedGraphToBigWig {input.bdg} {input.chromsize} {output}
  """


rule sortBed:
  input: "rlregions/counts.unsorted.bdg"
  output: "rlregions/counts.sorted.bdg"
  log: "logs/rlregions__sortBed.log"
  conda: "envs/bedgraph_to_bigwig.yaml"
  shell: """
    bedSort {input} {output}
  """


rule toBedGraph:
  input: "rlregions/counts.bed"
  output: "rlregions/counts.unsorted.bdg"
  shell: """
    awk 'BEGIN {{FS="\t"; OFS="\t"}} {{if ($5 > 0) print $1, $2, $3, $5}}' {input} > {output}
  """


rule bedtools_intersect:
  output: "rlregions/counts.bed"
  input:
    manifest = "rlregions/rlregion_manifest.csv",
    windows = "rlregions/10bp_windows.bed"
  log: "logs/rlregions__bedtools_intersect.log"
  conda: "envs/bedtools.yaml"
  shell: """
    intersectBed -a {input.windows} -b $(tr '\n' ' ' < {input.manifest}) -c | \
    awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, $2, $3, $4, $5, "."}}' - > {output}
  """


rule cat_peaks:
  """ This rule merges all the peaks in the peaksFinal/ folder into one file and sorts """
  input: "rlregions/rlregion_manifest.csv"
  output: "rlregions/cat_peaks.broadPeak"
  shell: """
    cat {input} | xargs cat | sort -k 1,1 -k2,2n > {output}
  """
  

rule get_rlregion_sample_manifest:
  input: "tables/rmap_samps.tsv"
  output: "rlregions/rlregion_manifest.csv"
  script: "scripts/buildRLRegionManifest.csv"


rule makewindows:
  input: "rlregions/hg38.chrom.sizes"
  output: "rlregions/10bp_windows.bed"
  log: "logs/rlregions__makewindows_10bp.log"
  conda: "envs/bedtools.yaml"
  shell: "bedtools makewindows -g {input} -w 10 -i srcwinnum | sort -k 1,1 -k2,2n > {output}"


rule download_chromsizes:
  output: "rlregions/hg38.chrom.sizes"
  shell: "wget -O {output} ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"


rule get_expression:
  input: quant,
  output: "tables/gene_expression.tsv"

rule get_rmap_samps:
  input: 
    manifest="rmap_manifest.csv",
    model="misc/fft_model.rda"
  output: "tables/rmap_samps.tsv"
  script: "scripts/compileRMapSamps.R"
  
rule build_model:
  input:
    manifest="rmap_manifest.csv",
    blacklist="misc/model_samples.csv",
    rlfs=rlfs_rdas
  output: "misc/fft_model.rda"
  script: "scripts/doBuildModel.R"
  
rule prep_model_samps:
  input:
    rmap="rmap_manifest.csv",
    rlfs=rlfs_rdas
  output: 
    mod_samps="misc/model_samples.csv"
  conda: "envs/jupyteR.yaml"
  notebook: "notebooks/filter_model_samples.r.ipynb"

rule get_pca:
  input: "tables/rmap_corr.tsv"
  output: "tables/rmap_pca.tsv"
  script: "scripts/doPCAAnalysis.R"

rule get_corr:
  input: coverage
  output: "tables/rmap_corr.tsv"
  script: "scripts/doCorrAnalysis.R"

rule gather_anno:
  input:  peaks
  output: "tables/rmap_anno.tsv"
  script: "scripts/doAnnoAnalysis.R"

rule prep_rlfs_out:
  """Aggregator for RLFS ahead of tarball"""
  input: peaks
  output: 
    tarball="misc/rmap_rlfs/tarball.txt",
    rdas=rlfs_rdas
  script: "scripts/doAnalyzeRLFS.R"

