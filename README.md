# RMapDB-Datasets

This `README` outlines all the steps required to fully build (or update) the datasets in RMapDB.
It has only been tested on `Ubuntu 18.04 LTS`. It will take several days/weeks to fully run the first
time, and it is not recommended without significant access to bioinformatics resources and a very high-speed
internet connection.

## Setup

1. Clone the repos

```shell
git clone https://github.com/Bishop-Laboratory/RMapDB-Datasets.git
git clone https://github.com/Bishop-Laboratory/RSeqCLI.git
git clone https://github.com/Bishop-Laboratory/RSeqR.git
```

2. Create the environment

```shell
cd RMapDB-Datasets/
conda install -c conda-forge mamba -y
mamba env create -f rmapdb-datasets.yml --force
conda activate rmapdb-datasets
```

3. Install RSeqCLI and RSeqR

```shell
pip install -e ../RSeqCLI/
R CMD build ../RSeqR
```

## Generate datasets

### Preliminary: Find Available genomes, RLFS bed files, etc

1. Download `QmRLFS-finder.py` (NOTE: First read `scripts/QmRLFS-finder/README.md`)

```shell
wget -O scripts/QmRLFS-finder/QmRLFS-finder.py https://raw.githubusercontent.com/piroonj/QmRLFS-finder/master/QmRLFS-finder.py
```

2. Generate the available genome information

```shell
Rscript scripts/makeAvailableGenomes.R
```

3. Generate RLFS Bed Files

```shell
Rscript scripts/makeRLFSBedFiles.R
```

4. (optional) Upload the results to AWS (Requires admin privileges)

```shell
aws s3 sync misc/ s3://rmapdb-data/misc/
aws s3 sync rlfs-beds/ s3://rmapdb-data/rlfs-beds/
```

### Sample peaks and coverage files

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

5. (optional) Upload datasets to aws (requires admin privledges)

```shell
aws s3 sync rmap-data/coverage/ s3://rmapdb-data/coverage/
aws s3 sync rmap-data/peaks/ s3://rmapdb-data/peaks/
aws s3 sync rmap-data/bam_stats/ s3://rmapdb-data/bam_stats/
aws s3 sync rmap-data/fq_stats/ s3://rmapdb-data/fq_stats/
```

6. (optional) Store .bam files (we used BOX for this part)

```shell
FTPSITE="ftp.box.com"
REMOTEDIR="RMapDB-Archive/bam/"
lftp -c "open -u user,pass $FTPSITE; mput -O $REMOTEDIR rmap-data/bam/" 
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

