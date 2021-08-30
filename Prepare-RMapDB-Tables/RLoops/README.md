# RLoop Table Analysis

The purpose of this analysis is to use the macs peaks generated via RSeq in order to determine a set of official R-loops which can be identified and named. This is necessary to establish an agreed-upon set for downstream analysis and comparison between studies. 

However, it is not yet clear what method for determining this RLoop set will be most appropriate. Three primary strategies have been devised:

(A) 
1. Generate 10 bp genomic windows
2. Construct a binary matrix of samples x windows (0 if sample doesn't have a peak in the window; 1 if sample does have a peak in the window)
3. Sum the matrix for every window to generate a bedGraph.
4. Call peaks on the bedGraph. 
5. Examine resulting peaks against RLFS

(B)
1. Generate 10 bp genomic windows
2. Construct a matrix of samples (i) X windows (j) in which cell i,j contains the -log10pval of sample i -- it is 1 if no overlap
3-5 same as (A) 3-5

(C) 
Same as (B) but using BPM instead from the .bw files. 

(D) 
1. Generate 10 bp genomic windows

