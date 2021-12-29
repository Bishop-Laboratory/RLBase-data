# RLBase-Datasets

This `README` outlines all the steps required to fully build (or update) the datasets in RLBase.
It has only been tested on `Ubuntu 20.04/18.04 LTS`. 
Processed data-sets can be accessed via the `RLHub` R package. 

**Time required**: It took around 3 weeks to fully run this protocol with 100 threads
and 250GiB of RAM allocated. Updating the database with new data will also take
substantial time and computing resources in proportion to the number and size 
of FASTQ files which must be aligned.

## Setup

1. Clone the repos

```shell
git clone https://github.com/Bishop-Laboratory/RLBase-data.git
git clone https://github.com/Bishop-Laboratory/RLPipes.git
git clone https://github.com/Bishop-Laboratory/RLSeq.git
```

2. Create the environment (requires [conda](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh))

```shell
cd RLBase-data/misc-data/
git clone https://github.com/srhartono/SkewR.git
wget https://github.com/KorfLab/StochHMM/archive/refs/tags/v0.38.tar.gz
tar -xvzf v0.38.tar.gz
cd StochHMM-0.38/
./configure
make
cd ../
RLBASEDIR=$(pwd)
export PATH="$RLBASEDIR/misc-data/StochHMM-0.38/:$PATH" 
conda install -c conda-forge mamba -y
mamba env create -f env.yml --force
conda activate rlbaseData
```

3. Install `RLPipes` and `RLSeq`

```shell
pip install -e ../RLPipes/
R -e "BiocManager::install(c('EnsDb.Hsapiens.v86', 'EnsDb.Mmusculus.v79'))"
R -e "install.packages(c('ggprism', 'tableone'), repos = 'http://cran.us.r-project.org')"
# If prior to bioc release:
# R -e "BiocManager::install(verion='devel')"
R -e "remotes::install_local('../RLHub/', dependencies=TRUE, force=TRUE)"
R -e "remotes::install_local('../RLSeq/', dependencies=TRUE, force=TRUE)"
```

## Generate datasets

The following steps taking to generate and update RLBase datasets. 
If you have write permissions on the AWS buckets used here, you can complete
every step if `awscli` is configured:

```shell
aws configure
```

If you do not, then you will only be able to mirror and download 
existing buckets.

### Preliminary

These preliminary steps illustrate how the original annotation data and 
genome info was collated.

1. Available genome info

This script downloads available genomes from the UCSC Genome Browser FTP.
It then uses `unique-kmers.py` (from the `khmer` pacakge) to calculate the 
effective genome size at various read lengths 
(36bp, 50bp, 75bp, 100bp, 125bp, 150p, 250bp). It then combines these results
with the other metadata available for the genome, saving the output in a
`data.frame` stored in `misc-data/available_genomes.rda`.

```shell
conda deactivate
Rscript scripts/makeAvailableGenomes.R
conda activate rlbaseData
```

2. RLFS Beds files

R-loop-forming sequences (RLFS) were obtained for each genome of interest using 
the `QmRLFS-finder.py` script. 

`QmRLFS-finder.py` was downloaded by the RLBase authors in October 2020 from the public repository [here](https://github.com/piroonj/QmRLFS-finder) and used with the permission of its creator. 

```shell
CORES=40  # Number of cores for parallel operations
Rscript scripts/makeRLFSBeds.R 1 FALSE  # Downloads needed genomes
Rscript scripts/makeRLFSBeds.R $CORES TRUE  # Runs QmRLFS-Finder.py in parallel
```

3. Generate Annotation Database

This script generates the genomic annotations needed by RLSuite. See the "Annotation Sources" section for further details.

<details>

  <summary>Annotation sources</summary>
  
Annotations relevant to R-loop biology were aggregated from a variety of sources:

  - UCSC tables (obtained via `rtracklayer::ucscTableQuery()`)
    - cpgIslandExt
    - centromeres
    - encodeCcreCombined
    - knownGene
    - microsat
    - rmsk
    - knownAl
    - wgRn
    - tRNAs 
    - coriellDelDup 
    - wgEncodeGencodePolyaV38
    - encRegTfbsClustered
  - UCSC Ensembl Gene GTFs
    - hg38 - [link](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz)
    - mm10 - [link](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ensGene.gtf.gz)
  - G4-Quadruplex Predictions - [link](https://figshare.com/ndownloader/files/6432597)
  - G4-Quadruplex Experiments - [link (+)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874_Na_K_PDS_plus_hits_intersect.bed.gz); [link (-)](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874_Na_K_PDS_minus_hits_intersect.bed.gz)
  - Encode 
    - Histone BED narrowPeak files (Full manifest in `misc-data/histone_encode_manifest.csv`)
    - RBP binding sites (ChIP and eCLiP) contributed by [this study](https://www.nature.com/articles/s41586-020-2077-3)
  - SkewR annotations of human and mouse genomes - generated by running `skewr` - [link](https://github.com/srhartono/SkewR)
  - RBP binding site predictions from oRNAment - [link](http://rnabiology.ircm.qc.ca/oRNAment/dashboard/)
  - Cohesin binding sites generated by the authors of [this study](https://doi.org/10.1093/nar/gkaa284).
  
</details>

```shell
# Runs the annotation script
CORES=15  # Set low due to high memory consumption per thread
Rscript scripts/getGenomicFeatures.R $CORES

# Compress resulting files
find misc-data/annotations/ -name "*.csv" -exec gzip -f {} \;
```

4. (optional) Upload the results to AWS (Requires admin privileges)

```shell
aws s3 cp misc-data/available_genomes.rda s3://rlbase-data/misc/
aws s3 cp misc-data/rlfs/ s3://rlbase-data/rlfs-beds/ --recursive --exclude "*" --include "*.bed" --exclude "*/*"
aws s3 sync misc-data/annotations/ s3://rlbase-data/annotations/
```

## Run RLPipes on all public samples

This part of the protocol requires a catalog of publicly-available R-loop-mapping
samples, hand-curated in the manner described in the RLSuite publication. 

The following steps are performed to generate (or update) the database:

1. Prepare catalog of publicly-available samples. The [current catalog](https://github.com/Bishop-Laboratory/RLBase-data/raw/main/rlbase-data/rlbase_catalog.xlsx) can serve as a template.

2. Make pipeline manifests from catalog

```shell
CATALOG="rlbase-data/rlbase_catalog.xlsx"
MANIFEST="rlbase-data/rlbase_manifest.csv"
Rscript scripts/makeManifest.R $CATALOG $MANIFEST
```

3. Build the pipeline config
```shell
RLPIPESOUT="rlbase-data/rlpipes-out/"
RLPipes build $RLPIPESOUT $MANIFEST
cp $RLPIPESOUT/config.tsv $RLPIPESOUT/config.save.tsv
```

4. Wrangle the config (fix genomes and add condition types)

```shell
Rscript scripts/wrangleConfig.R $CATALOG $RLPIPESOUT/config.tsv
```

5. Remove any datasets which have to be re-run

```shell
Rscript scripts/cleanOld.R $CATALOG $RLPIPESOUT/config.tsv
```

6. Then check the pipeline to ensure it will work correctly:

```shell
RLPipes check $RLPIPESOUT --bwamem2 --noreport --tsv
```

7. Run the pipeline

```shell
CORES=44
RLPipes run $RLPIPESOUT --bwamem2 --noreport --tsv -t $CORES
```

8. Upload datasets to aws (requires admin privledges)

```shell
cd $RLBASEDIR
aws s3 sync --size-only $RLPIPESOUT/coverage/ s3://rlbase-data/coverage/
aws s3 sync --size-only $RLPIPESOUT/peaks_save/ s3://rlbase-data/peaks/
aws s3 sync --size-only $RLPIPESOUT/quant_save/ s3://rlbase-data/quant/
aws s3 sync --size-only $RLPIPESOUT/bam_stats/ s3://rlbase-data/bam_stats/
aws s3 sync --size-only $RLPIPESOUT/fastq_stats/ s3://rlbase-data/fastq_stats/
```

9. Store .bam files (we used BOX for this part)

```shell
FTPSITE="ftp.box.com"
REMOTEDIR="RLBase-Archive/"
lftp -c "open -u $(read -p "User: ";echo $REPLY),$(read -sp "Password: ";echo $REPLY) $FTPSITE; mirror -P 4 -n -R $RLPIPESOUT/bam/ $REMOTEDIR"
```

10. Calculate RLFS enrichment for each sample

```shell
CORES=44
Rscript scripts/rlfsAnalyze.R $RLPIPESOUT/peaks $RLPIPESOUT/rlfs_rda $CORES
```

11. Upload to AWS

```shell
aws s3 sync --size-only $RLPIPESOUT/rlfs_rda/ s3://rlbase-data/rlfs_rda/
```

## Downstream data processing

0. If you didn't do the previous steps, then download the data needed for this part:

```shell
RLPIPESOUT="rlbase-data/rlpipes-out/"
aws s3 sync s3://rlbase-data/peaks/ $RLPIPESOUT/peaks/ 
aws s3 sync s3://rlbase-data/rlfs_rda/ $RLPIPESOUT/rlfs_rda/ 
aws s3 sync s3://rlbase-data/quant/ $RLPIPESOUT/quant/ 
cd $RLPIPESOUT/quant/ && find . -name "*.tar.xz" -exec tar -xJf {} \; && cd $RLBASEDIR
aws s3 sync s3://rlbase-data/coverage/ $RLPIPESOUT/coverage/ 
aws s3 sync s3://rlbase-data/fastq_stats/ $RLPIPESOUT/fastq_stats/ 
aws s3 sync s3://rlbase-data/bam_stats/ $RLPIPESOUT/bam_stats/ 
```

### Build discriminator model

1. Run the interactive model-building app (requires AWS privileges)

```shell
CONFIG=$RLPIPESOUT/config.tsv
HOST="0.0.0.0"
PORT=4848
Rscript scripts/models.R $CONFIG $HOST $PORT
```

2. Classify all samples

```shell
Rscript scripts/classifySamples.R
```

### Get R-loop consensus

1. Select final peakset

```shell
Rscript scripts/prepareConsensus.R
```

2. Run the R-loop consensus pipeline

```shell
CORES=44
MANIFEST="rlbase_samples.tsv"
BLACKLIST="https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
snakemake --snakefile scripts/rlregions.smk -d rlbase-data/ --cores $CORES --dryrun --config manifest=$MANIFEST blacklist=$BLACKLIST # Verify
snakemake --snakefile scripts/rlregions.smk -d rlbase-data/ --cores $CORES --config manifest=$MANIFEST blacklist=$BLACKLIST --dag | dot -Tpng > rlbase-data/rlregions/dag.png  # DAG
snakemake --snakefile scripts/rlregions.smk -d rlbase-data/ --config manifest=$MANIFEST blacklist=$BLACKLIST --cores $CORES # Run it
```

### Update processed datasets

1. Annotate peaks (genomic features and genes)

```shell
CORES=44
PEAKS=$RLPIPESOUT/peaks
OUTANNO="rlbase-data/misc/annotatedPeaks.tsv"
Rscript scripts/annotatePeaks.R $PEAKS $OUTANNO $CORES
```

2. Build the new `gene_expression` table and calculate correlations. 

```shell
# Creates misc/gene_expression.rda
MANIFEST="rlbase-data/rlbase_samples.tsv"
GENE_EXP_TABLE="rlbase-data/misc/gene_expression.rda"
QUANTDIR="rlbase-data/rlpipes-out/quant/"
Rscript scripts/buildExpression.R $QUANTDIR $GENE_EXP_TABLE $MANIFEST
```

3. Get the genomic feature -> rlregion mapping

```shell
# Creates misc/rlregion_features.rda
Rscript scripts/rlregionsToFeatures.R
```

4. Build the new correlation matrix. 

```shell
Rscript scripts/gsgCorr.R
```

5. Get RLRegion count matrix

```shell
CORES=44
Rscript scripts/rlregionCountMat.R $CORES 1
```

6. Rebuild all RLRanges and HTML notebooks

```shell
CORES=32
Rscript scripts/runRLSeq.R $CORES
```

7. Upload results to AWS
```shell
aws s3 sync misc-data/reports s3://rlbase-data/reports
aws s3 sync misc-data/rlranges s3://rlbase-data/rlranges
```


8. Build/Update the RLHub 

```shell
Rscript scripts/prepRLHub.R
find misc-data/rlhub/ -name "*.rda" -exec aws s3 cp {} s3://rlbase-data/RLHub/ \;
```

9. Update the RLHub Genome Browser TrackHub

```shell
Rscript scripts/buildGenomeBrowserHub.R
aws s3 sync misc-data/RLBase_TrackHub/ s3://rlbase-data/misc/RLBase_TrackHub/
```


