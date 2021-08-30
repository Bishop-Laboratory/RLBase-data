import pandas as pd
from os import path
import math

## Constants ##
RUN_SHEET = "run_sheet.csv"
debug = False
if debug:
  RUN_SHEET = "run_sheet.small.csv"
###############

## Read config ##
config = pd.read_csv(RUN_SHEET)
sample=config['sample']
experiments=[exp.split(",") for exp in config['experiment']]
genome=config['genome']
paired_end=config['paired_end']
#################

## Helpers ##
def test_pe(wildcards):
  return [paired_end[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]

def find_sra_replicates(wildcards):
  srr_list = [experiments[idx] for idx, element in enumerate(sample) if element == wildcards.sample][0]
  replicates = expand("fastqs_raw/{sample}/{srr_acc}.fastq",
         srr_acc=srr_list, sample=wildcards.sample)
  return replicates
  
  
def pe_test_fastp(wildcards):
  pe = test_pe(wildcards)
  if pe:
      res="--interleaved_in "
  else:
      res=""
  return res


def get_fq_salmon(wildcards):
  pe = test_pe(wildcards)
  if pe:
      res="-1 fastqs_prepped/" + wildcards.sample + "." + wildcards.genome +\
      ".R1.fastq -2 fastqs_prepped/" + wildcards.sample + "." +\
      wildcards.genome + ".R2.fastq"
  else:
      res="-r fastqs_prepped/" + wildcards.sample + "." + wildcards.genome + ".R1.fastq"
  return res

  
############


## Rules ##

rule output:
  input: expand("quant/{sample}_{genome}/quant.sf", zip, sample=sample, genome=genome)


rule salmon:
  input: 
    fq="fastqs_prepped/{sample}.{genome}.R1.fastq",
    ind="salmon_ind/{genome}/versionInfo.json"
  output: "quant/{sample}_{genome}/quant.sf"
  params:
    pe_fq=get_fq_salmon,
    ind="salmon_ind/{genome}",
    outdir="quant/{sample}_{genome}"
  log: "logs/quant/{sample}_{genome}__quant.log"
  conda: "envs/salmon.yaml"
  threads: 6
  priority: 11
  shell:
    """
    (
      salmon quant -i {params.ind} -l A {params.pe_fq} --validateMappings -o {params.outdir} -p {threads}
    ) &> {log}
    
    """



rule download_annotations:
  output: 
    fasta="salmon_ind/hg38/versionInfo.json"
  params:
    name="salmon_ind/hg38"
  shell: """
    wget -O default.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/salmon_sa_index?tag=default
    tar -xvzf default.tgz
    rm default.tgz
    rm -rf {params.name}
    mv default {params.name}
  """  


rule cleanup_fq:
  input: "fastqs_trimmed/{sample}.{genome}__trimmed.fastq"
  output: "fastqs_prepped/{sample}.{genome}.R1.fastq"
  conda: "envs/sratools.yaml"
  log: "logs/cleanup_fq/{sample}.{genome}__cleanup_fq.log"
  params:
    ispe=test_pe,
    outdir="fastqs_prepped/"
  shell: """
  (
    if [ {params.ispe} == "True" ]; then
      reformat.sh ow=t in={input} out1={params.outdir}/{wildcards.sample}.{wildcards.genome}.R1.fastq \
      out2={params.outdir}/{wildcards.sample}.{wildcards.genome}.R2.fastq
    else
      mv {input} {output}
    fi
  ) &> {log}
  """
  

# TODO: Pipe this part
rule fastp:
  input: "fastqs_merged/{sample}.{genome}__merged.fastq"
  output:
      trimmed=temp("fastqs_trimmed/{sample}.{genome}__trimmed.fastq"),
      html="QC/fastq/html/{sample}.{genome}.html",
      json="QC/fastq/json/{sample}.{genome}.json"
  conda: "envs/fastp.yaml"
  log: "logs/fastp/{sample}.{genome}__fastp_pe.log"
  priority: 10
  params:
      extra=pe_test_fastp
  threads: 4
  shell: """
  (fastp -i {input} --stdout {params.extra}-w {threads} -h {output.html} -j {output.json} > {output} ) &> {log}
  """

rule merge_replicate_reads:
  input: find_sra_replicates
  output: temp("fastqs_merged/{sample}.{genome}__merged.fastq")
  log: "logs/merge_fastq/{sample}_{genome}__merge_fastq.log"
  shell: """
  (cat {input} > {output}) &> {log}
  """
  
  
# If debugging, will use a less stable but much faster version of this rule which outputs a small fastq file
# This should always produced interleaved fq even if paired end. I don't like this solution, but it should hold.
# reformat.sh (from BBTools) will interleave the paired-end files
rule sra_to_fastq:
  input: "sras/{sample}/{srr_acc}/{srr_acc}.sra"
  output: temp("fastqs_raw/{sample}/{srr_acc}.fastq")
  conda: "envs/sratools.yaml"
  threads: 1
  log: "logs/sra_to_fastq/{sample}_{srr_acc}__sra_to_fastq_pe.log"
  params:
      output_directory="sras/{sample}/",
      fqdump="--skip-technical --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --split-3 "
  shell: """(
  cd {params.output_directory}
  fastq-dump {params.fqdump}-O ../../fastqs_raw/{wildcards.sample}/ {wildcards.srr_acc}
  cd ../../fastqs_raw/{wildcards.sample}/
  if test -f {wildcards.srr_acc}_2.fastq; then
      echo "Paired end -- interleaving"
      reformat.sh in1={wildcards.srr_acc}_1.fastq in2={wildcards.srr_acc}_2.fastq out={wildcards.srr_acc}.fastq overwrite=true
      rm {wildcards.srr_acc}_1.fastq && rm {wildcards.srr_acc}_2.fastq
  else
      echo "Single end -- finished!"
  fi
  ) &> {log}
  """


# TODO: Figure out how to use the pipes
# TODO: Probably need to specify this version of prefetch and/or find alternative to it...
# TODO: Retry if fails due to network error
rule download_sra:
  output: "sras/{sample}/{srr_acc}/{srr_acc}.sra"
  conda: "envs/sratools.yaml"
  log: "logs/download_sra/{sample}__{srr_acc}__download_sra.log"
  params:
      output_directory = "sras/{sample}/"
  threads: 10
  shell: """
          (
          cd {params.output_directory}
          prefetch {wildcards.srr_acc} -f yes
          ) &> {log}
          """


            
