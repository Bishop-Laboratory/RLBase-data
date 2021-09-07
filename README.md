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
cd RLBase-data/
RLBASEDIR=$(pwd)
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

4. Fix any genomes which are mislabeled

```shell
Rscript scripts/fixGenomes.R $CATALOG $RLPIPESOUT/config.tsv
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

## Build discriminator model

1. Calculate RLFS enrichment for each sample

```shell
CORES=44
Rscript scripts/rlfsAnalyze.R $RLPIPESOUT/peaks $RLPIPESOUT/rlfs_rda $CORES
```

2. Upload to AWS

```shell
aws s3 sync $RLPIPESOUT/rlfs_rda/ s3://rlbase-data/rlfs_rda/
```

3. Get the samples for model buidling

```shell
CONFIG="rlbase-data/rlpipes-out/config.tsv"
HOST="0.0.0.0"
PORT=4848
Rscript scripts/selectSamples.R $CONFIG $HOST $PORT
```

4. Then build the model

```shell
Rscript scripts/buildModel.R
```

5. Re-build RLSeq with new data

```shell
cp misc-data/FFTModel.rda ../RLSeq/data/
cp misc-data/prepFeatures.rda ../RLSeq/data/
R -e "remotes::install_local('../RLSeq/', dependencies=TRUE, force=TRUE)"
```

6. Classify samples

```shell
Rscript script/classifySamples.R
```


### Sample peaks and coverage files

```shell
snakemake --snakefile rlbase-datasets.smk --configfile ../RSeqCLI/tests/rseq_out_public_rna/config.json -d rmap-data/ --notebook-listen 0.0.0.0:6123 --edit-notebook misc/rmap_blacklist.csv
```

1. Create the RMap Sample manifest from the latest hand-made excel manifest

```shell
Rscript scripts/createRMapDBManifests.R
```

2. Build the RSeqCLI config file for R-loop mapping samples (this will take several minutes)

```shell
RSeqCLI build rmap-data/rmap/ db-data/rmap_manifest.csv
```

3. Check that the workflow will execute properly

```shell
# --no-report indicates report files will not be generated
RSeqCLI check rmap-data/rmap/ --bwamem2 -t 44 --no-report
```

4. Execute the workflow (will take several days/weeks to complete from scratch)

```shell
# --no-report indicates report files will not be generated
RSeqCLI run rmap-data/rmap/ --bwamem2 -t 44 --no-report
```



### Rebuild datasets for RSeqR

1. Run the RLFS analysis tool from RSeqR on all data

```shell
Rscript scripts/runRLFSAnalysis.R
```

2. Run the model-generation script (credit to Daniel Montemayor)

```shell
# Creates misc/rdata/model.rda
Rscript buildModel.R
```

3. Use the new model/manifest to build the new `rmapdb_samples` table. 

```shell
# Creates misc/rdata/rmapdb_samples.rda
Rscript buildRMapDBSamples.R  
```

4. Build the new RMapDB table and peaks to build the `rmapdb_genfeat` table.

```shell
# Creates misc/rdata/rmapdb_genfeat.rda
Rscript buildRMapDBGenFeat.R  
```

5. Use new RMapDB samps table to build the new `rloop_regions` table

```shell
# Creates misc/rdata/rloop_regions.rda
Rscript buidRLoopRegiones.R
```

6. Use preceding to build the new `rloop_signal` table

```shell
# Creates misc/rdata/rloop_signal.rda
Rscript buildRloopSignal.R
```

7. Build the new `gene_expression` table. 

```shell
# Creates misc/rdata/gene_expression.rda
Rscript buildGeneExpression.R
```

8. Update the data in the RSeqR package and rebuild RSeqR

```shell
cp -rf misc/rdata ../RSeqR/data/
R CMD build ../RSeqR
```

### Generate reports, upload, and clean up

1. Rerun RSeqCLI to generate the reports using the latest RSeqR version

```shell
RSeqCLI run rmap-data/rmap/ --bwamem2 -t 44
```

2. Generate the new R-loops bigWig/bed file

```shell
# Creates misc/rloop_regions.bw & misc/rloop_regions.bed
Rscript buildRLRegions.bw.R
```

3. Generate data dump

```shell
# Creates misc/rmapdb_peaks.tar.xz & misc/rmapdb_coverage.tar.xz & misc/rmapdb_reports.tar.xz
Rscript makeDataDumps.R
```

4. Upload data to AWS 

```shell
aws s3 sync rmap-data/report/data/ s3://rmapdb-data/report/data/
aws s3 sync rmap-data/report/html/ s3://rmapdb-data/report/html/
aws s3 sync misc/ s3://rmapdb-data/misc/
```

