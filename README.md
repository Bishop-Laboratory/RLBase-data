# RMapDB-Datasets

This `README` outlines all the steps required to fully build (or update) the datasets in RMapDB.
It has only been tested on `Ubuntu 18.04 LTS`. It will take several days/weeks to fully run the first
time, and it is not recommended without access to high-performance computing resources.

## Setup

1. Clone the repos

```shell
git clone https://github.com/Bishop-Laboratory/RLBase-data.git
git clone https://github.com/Bishop-Laboratory/RLPipes.git
git clone https://github.com/Bishop-Laboratory/RLSeq.git
```

2. Create the environment

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

3. Install RLPipes and RLSeq

```shell
pip install -e ../RLPipes/
R -e "BiocManager::install(c('EnsDb.Hsapiens.v86', 'EnsDb.Mmusculus.v79'))"
R -e "install.packages('ggprism', repos = 'http://cran.us.r-project.org')"
R -e "remotes::install_local('../RLSeq/', dependencies=TRUE, force=TRUE)"
```

## Generate datasets

The following details the steps taking to generate and update RLBase datasets. 
If you have write permissions on the AWS buckets used here, you can complete
every step if `awscli` is configured:

```shell
aws configure
```

If you do not, then you will only be able to mirror and download 
existing buckets.

### Preliminary

These preliminary steps just illustrate how the original annotation data and 
genome info was collated.

#### Available genome info

Generate the available genome information using the script below (takes several hours):

```shell
conda deactivate
R -e "source('scripts/makeAvailableGenomes.R')"
conda activate rlbaseData
```

#### RLFS Beds files

These data are now permanently available from the bucket located
here: s3://rlbase-data/rlfs-beds/

To download all files, run the following command:

```shell
aws s3 sync s3://rlbase-data/rlfs-beds/ rlfs-beds/
```

These files were generated via the following:

1. Download `QmRLFS-finder.py` (NOTE: First read `scripts/QmRLFS-finder/README.md`)

```shell
wget -O scripts/QmRLFS-finder/QmRLFS-finder.py https://raw.githubusercontent.com/piroonj/QmRLFS-finder/master/QmRLFS-finder.py
```

2. Download and unpack genomes 

```shell
Rscript scripts/makeRLFSBeds.R 1 FALSE
```

3. Generate RLFS Bed Files (Takes several hours, depending on parallelization and download speed)

```shell
CORES=40  # Number of cores for parallel operations
Rscript scripts/makeRLFSBeds.R $CORES TRUE
```

4. Generate SkewR tracks

```shell
cd misc-data/SkewR/
perl bin/RunGC-SKEW.pl -s ~/.rseq_genomes/hg38/hg38.fa -m model/GC_SKEW_7600.hmm -g ~/.rseq_genomes/hg38/hg38.ensGene.bed -b ~/.rseq_genomes/hg38/hg38.cpg.bed -o skewr_res_hg38 -z 20
```

4. (optional) Upload the results to AWS (Requires admin privileges)

```shell
aws s3 sync misc-data/ s3://rmapdb-data/misc/
aws s3 sync rlfs-beds/ s3://rmapdb-data/rlfs-beds/
```

## Run RLPipes on all public samples

This part of the protocol requires a catalog of publicly-available R-loop-mapping
samples, hand-curated in the manner described in the RLSuite publication. 

The previous catalog is available [here](https://github.com/Bishop-Laboratory/RLBase-data/raw/main/rlbase-data/rlbase_catalog.xlsx)

The following are the steps performed to generate the database:

1. Prepare catalog of publicly-available samples. Current catalog can serve as a template.

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
CORES=179
RLPipes run $RLPIPESOUT --bwamem2 --noreport --tsv -t $CORES
```

8. Archive quant folders

```shell
cp -r $RLPIPESOUT/quant $RLPIPESOUT/quant_save
cd $RLPIPESOUT/quant_save
find . -type d -maxdepth 1 -mindepth 1 -exec tar cfJ {}.tar.xz {} \;
find . -type d -maxdepth 1 -mindepth 1 -exec rm -rf {} \;
cd $RLBASEDIR
```

9. Archive peaks

```shell
mkdir $RLPIPESOUT/peaks_save
find $RLPIPESOUT/peaks -type f -maxdepth 1 -mindepth 1 -name "*.broadPeak" -exec cp {} $RLPIPESOUT/peaks_save \;
```

10. Make tarballs

```shell
cd $RLPIPESOUT
tar cfJ peaks.tar.xz peaks/
tar cfJ coverage.tar.xz coverage/
tar cfJ quant.tar.xz quant/
tar cfJ logs.tar.xz logs/
tar cfJ fastq_stats.tar.xz fastq_stats/
tar cfJ bam_stats.tar.xz bam_stats/
mkdir datadump
find . -type f -maxdepth 1 -name "*.tar.xz" -exec mv {} datadump/ \;
cp config.tsv datadump/
cp config.json datadump/
cp ../rlbase_catalog.xlsx datadump/
cp ../rlbase_manifest.csv datadump/
```

11. Upload datasets to aws (requires admin privledges)

```shell
aws s3 sync $RLPIPESOUT/coverage/ s3://rlbase-data/coverage/
aws s3 sync $RLPIPESOUT/peaks_save/ s3://rlbase-data/peaks/
aws s3 sync $RLPIPESOUT/quant_save/ s3://rlbase-data/quant/
aws s3 sync $RLPIPESOUT/bam_stats/ s3://rlbase-data/bam_stats/
aws s3 sync $RLPIPESOUT/fastq_stats/ s3://rlbase-data/fastq_stats/
aws s3 sync $RLPIPESOUT/datadump/ s3://rlbase-data/datadump/
```

12. Store .bam files (we used BOX for this part)

```shell
FTPSITE="ftp.box.com"
REMOTEDIR="RMapDB-Archive/"
lftp -c "open -u $(read -p "User: ";echo $REPLY),$(read -sp "Password: ";echo $REPLY) $FTPSITE; mirror -P 4 -n -R $RLPIPESOUT/bam/ $REMOTEDIR"
```

13. Calculate RLFS enrichment for each sample

```shell
CORES=44
Rscript scripts/rlfsAnalyze.R $RLPIPESOUT/peaks $RLPIPESOUT/rlfs_rda $CORES
```

14. Upload to AWS

```shell
aws s3 sync $RLPIPESOUT/rlfs_rda/ s3://rlbase-data/rlfs_rda/
```

## Build discriminator model

0. If you didn't do the above, then download the data needed for this:

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


1. Get the samples for model buidling

TODO: NEED TO Add condition key for NEG, POS, NULL
TODO: Need to highlight the new samples and have another column previously determined vs new
TODO: Need a reference for judging whether to discard -- show the ideals
TODO: Maybe use a PCA instead of the fourier transform graph
TODO: Maybe take out the annotation filter

```shell
CONFIG=$RLPIPESOUT/config.tsv
HOST="0.0.0.0"
PORT=4848
Rscript scripts/selectSamples.R $CONFIG $HOST $PORT
```

2. Then build the model

```shell
Rscript scripts/buildModel.R
```

3. Re-build RLSeq with new data

```shell
cp misc-data/fftModel.rda ../RLSeq/data/
cp misc-data/prepFeatures.rda ../RLSeq/data/
R -e "remotes::install_local('../RLSeq/', dependencies=TRUE, force=TRUE)"
```

4. Classify samples

```shell
Rscript scripts/classifySamples.R
```

## Get R-loop consensus

1. Select final peakset

```shell
Rscript scripts/prepareConsensus.R
```

3. Run the R-loop consensus pipeline

```shell
CORES=44
MANIFEST="rlbase_samples.tsv"
BLACKLIST="https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
snakemake --snakefile scripts/rlregions.smk -d rlbase-data/ --cores $CORES --dryrun --config manifest=$MANIFEST blacklist=$BLACKLIST # Verify
snakemake --snakefile scripts/rlregions.smk -d rlbase-data/ --cores $CORES --config manifest=$MANIFEST blacklist=$BLACKLIST --dag | dot -Tpng > misc-data/rlregions.png  # DAG
snakemake --snakefile scripts/rlregions.smk -d rlbase-data/ --config manifest=$MANIFEST blacklist=$BLACKLIST --cores $CORES # Run it
```

## Other

0. Download the cohesin peaks

```shell
aws s3 sync s3://rlbase-data/misc/cohesin_peaks/ misc-data/cohesin_peaks/ 
```

1. Compile genomic features

```shell
CORES=15  # Set low due to high memory consumption per thread
Rscript scripts/getGenomicFeatures.R $CORES
```

2. Upload to AWS (Optional)

```shell
aws s3 sync misc-data/annotations/ s3://rlbase-data/misc/annotations/ 
```

2. Annotate peaks (genomic features and genes)

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
GENE_EXP_TABLE="rlbase-data/misc/gene_expression.csv"
QUANTDIR="rlbase-data/rlpipes-out/quant/"
Rscript scripts/buildExpression.R $QUANTDIR $GENE_EXP_TABLE $MANIFEST
```

3. Get the genomic feature -> rlregion mapping

```shell
# Creates misc/rlregion_features.rda
Rscript scripts/rlregionsToFeatures.R
```

3. Build the new correlation matrix. 

4. Build tables for DB

5. Update RLSeq with new data

6. Generate reports



4. Build the SQL database (optional)
























