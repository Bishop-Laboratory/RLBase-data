## Get the RLFS download and fix the width/score issue that keeps UCSC from displaying them correctly
RLFS_PEAKS <- "hg38.rlfs.bed"
RLFS_PEAKS_URI <- "s3://rmapdb-data/misc/rlfs.hg38.fixed.bed"
download.file("https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/hg38.rlfs.bed", 
              destfile = "misc/hg38.rlfs.bed")
rlfs <- rtracklayer::import("misc/hg38.rlfs.bed")
names(rlfs) <- paste0("hg38_", seq(width(rlfs)))
score(rlfs) <- 0
rtracklayer::export(rlfs, con = RLFS_PEAKS)
system(paste0("aws s3 --region us-west-2 cp ", RLFS_PEAKS, " ", RLFS_PEAKS_URI, " --acl public-read"))