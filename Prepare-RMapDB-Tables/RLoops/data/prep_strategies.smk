import os
import glob
import re
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

################################# MAGIC VALUES #################################
winsize = 10  # Set lower for higher accuracy but longer runtimes
S3_BW_URI = "rmapdb-data/bigwigs/rseq-coverage-unstranded/"
DEEPTOOLS_THREADS = 44
PEAKS_HG38_DIR = "peaksHG38/"
PEAKS_ANALYZE_DIR = "peaksFinal/"
################################################################################

def find_peaks(wildcards):
  files = glob.glob(PEAKS_ANALYZE_DIR + '*.broadPeak', recursive=True)  
  return files

def find_peaks2(wildcards):
  files = glob.glob(PEAKS_HG38_DIR + '*.broadPeak', recursive=True)  
  return files

def get_rl_signal(wildcards):
  files = find_peaks2(wildcards)
  infiles = [ 
    S3_BW_URI +\
    re.sub("^" + PEAKS_HG38_DIR + "([A-Za-z0-9_\-]+)_hg38.unstranded.broadPeak$", 
           "\\1", file) +\
    ".hg38.bw" for file in files 
  ]
  # Problematic file... Removing.
  infiles.remove("rmapdb-data/bigwigs/rseq-coverage-unstranded/peaksHG38/SRX2675011_HKE293-ΔHC-V5ChIP_hg38.unstranded.broadPeak.hg38.bw")
  return(
    # From www.snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html
    S3.remote(infiles, keep_local=True)
  )
  

################################# SNAKE RULES ##################################

rule output:
  input:
    strategy_a_peak="rloops/strategy.a__" + str(winsize) + "bp__peaks.narrowPeak",
    strategy_a_counttsv="rloops/strategy.a__" + str(winsize) + "bp__peaks.countMat.tsv.xz",
    strategy_a_bw="counts/strategy.a__" + str(winsize) + "bp__counts.bw",
    # strategy_b_peak="rloops/strategy.b__" + str(winsize) + "bp__peaks.narrowPeak",
    # strategy_b="counts/strategy.b__" + str(winsize) + "bp__counts.bw",
    # strategy_c_peak="rloops/strategy.c__" + str(winsize) + "bp__peaks.narrowPeak",
    # strategy_c="counts/strategy.c__" + str(winsize) + "bp__counts.bw",
    rloop_signal="rloop_signal.tsv.xz"
    
    
    
rule bedtools_countTSV_stratA:
  output: "rloops/{strat}__{winsize}bp__peaks.countMat.tsv"
  input:
    rloops = "rloops/{strat}__{winsize}bp__peaks.narrowPeak"
  log: "logs/{strat}__{winsize}bp__bedtools_intersect.log"
  conda: "envs/bedtools.yaml"
  params:
    peaks_analyze=PEAKS_HG38_DIR
  shell: """
    intersectBed -a {input.rloops} -b {params.peaks_analyze}/* -wao | \
    awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print $4, $15, $18, $20, "."}}' - > {output}
  """
    
    
rule rloop_signal:
  input:
    signals=get_rl_signal,
    peaks='rloops.bed.xz'
  output:
    npz="rloop_signal.npz",
    tsv="rloop_signal.tsv"
  conda: "envs/deeptools.yaml"
  threads: DEEPTOOLS_THREADS
  log: "logs/deeptools.log"
  params:
    bed="rloops.bed"
  shell: """
  (
    unxz {input.peaks}
    multiBigwigSummary BED-file -b {input.signals} --BED {params.bed} \
    -p {threads} -v -o {output.npz} --outRawCounts {output.tsv}
  ) &> {log}
  """
    
    
rule xz_file:
  input: "{filename}"
  output: "{filename}.xz"
  shell: "xz {input}"


# TODO: Refine peak

rule callpeak_from_bdg:
  input: "counts/{strat}__{winsize}bp__counts.sorted.bdg"
  output: "rloops/{strat}__{winsize}bp__peaks.narrowPeak"
  log: "logs/{strat}__{winsize}bp__callpeak_from_bdg.log"
  conda: "envs/macs.yaml"
  shell: """
    macs3 bdgpeakcall -i {input} -o {output}
  """

# TODO: Constrain by blacklist


rule toBigWig:
  input:
    bdg="counts/{strat}__{winsize}bp__counts.sorted.bdg",
    chromsize="hg38.chrom.sizes"
  output: "counts/{strat}__{winsize}bp__counts.bw"
  log: "logs/{strat}__{winsize}bp__toBigWig.log"
  conda: "envs/bedgraph_to_bigwig.yaml"
  shell: """
    bedGraphToBigWig {input.bdg} {input.chromsize} {output}
  """


rule sortBed:
  input: "counts/{strat}__{winsize}bp__counts.unsorted.bdg"
  output: "counts/{strat}__{winsize}bp__counts.sorted.bdg"
  log: "logs/{strat}__{winsize}bp__sortBed.log"
  conda: "envs/bedgraph_to_bigwig.yaml"
  shell: """
    bedSort {input} {output}
  """


rule toBedGraph:
  input: "counts/{strat}__{winsize}bp__counts.bed"
  output: "counts/{strat}__{winsize}bp__counts.unsorted.bdg"
  shell: """
    awk 'BEGIN {{FS="\t"; OFS="\t"}} {{if ($5 > 0) print $1, $2, $3, $5}}' {input} > {output}
  """


rule bedtools_intersect_stratA:
  output: "counts/strategy.a__{winsize}bp__counts.bed"
  input:
    windows = "windows/{winsize}bp_windows.bed"
  log: "logs/strategy.a__{winsize}bp__bedtools_intersect.log"
  conda: "envs/bedtools.yaml"
  params:
    peaks_analyze=PEAKS_ANALYZE_DIR
  shell: """
    intersectBed -a {input.windows} -b {params.peaks_analyze}/* -c | \
    awk 'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, $2, $3, $4, $5, "."}}' - > {output}
  """


rule bedtools_map_stratB:
  input:
    windows="windows/{winsize}bp_windows.bed",
    peakcat="cat_peaks/cat_peaks.broadPeak"
  log: "logs/strategy.b__{winsize}bp__mapBed.log"
  output: "counts/strategy.b__{winsize}bp__counts.bed"
  conda: "envs/bedtools.yaml"
  shell: """
    mapBed -c 7 -o mean -b {input.peakcat} -a {input.windows} -null 0 > {output}
  """


rule bedtools_map_stratC:
  input:
    windows="windows/{winsize}bp_windows.bed",
    peakcat="cat_peaks/cat_peaks.broadPeak"
  log: "logs/strategy.c__{winsize}bp__mapBed.log"
  output: "counts/strategy.c__{winsize}bp__counts.bed"
  conda: "envs/bedtools.yaml"
  shell: """
    mapBed -c 9 -o mean -b {input.peakcat} -a {input.windows} -null 0 > {output}
  """


rule cat_peaks:
  """ This rule merges all the peaks in the peaksFinal/ folder into one file and sorts """
  output: "cat_peaks/cat_peaks.broadPeak"
  params:
    peaks_analyze=PEAKS_ANALYZE_DIR
  shell: """
    find {params.peaks_analyze}/ -name "*.broadPeak" | xargs cat | sort -k 1,1 -k2,2n > {output}
  """


rule makewindows:
  input: "hg38.chrom.sizes"
  output: "windows/{winsize}bp_windows.bed"
  log: "logs/makewindows__{winsize}bp.log"
  conda: "envs/bedtools.yaml"
  shell: "bedtools makewindows -g {input} -w {wildcards.winsize} -i srcwinnum | sort -k 1,1 -k2,2n > {output}"


rule download_chromsizes:
  output: "hg38.chrom.sizes"
  shell: "wget -O {output} ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"

