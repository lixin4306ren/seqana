# seqana
chip-seq, bs-seq, rna-seq processing scripts

## installation  
```
install.packages("devtools")
library("devtools")
install_github("lixin4306ren/seqana")
```

## chip-seq related functions
1. `diff_peaks_replicate`   
merge peaks based on replicates and find differentially peaks between 2 samples
2. `merge_peaks`  
merge neighbouring peaks to large one by given distance
3. `get_common_peak3`  
extract common peaks based on replicates
4. `cal.rpkm.chip`  
this function calculate RPKM for chip-seq data, support bam and bed input format
5. `annotatePeaks`  
find nearest gene asscosiated with given peaks, deault database is human.gencodeV19

## bs-seq related functions
1. `plot_meth_sliding_bsseq`  
plot and calculate methylation along sliding windows from bsseq object
2. `meth_count_region_for_bsseq`  
calculate methylation level based on given Grange from bsseq object

## rna-seq related funtions
## utility functions
1. `generate.sliding.window`  
generate sliding window ranges
2. `bedtools.coveragebed.range`    
calculate overlap % between 2 bed files, need bedtools installed
3. `bedtools.coveragebed`  
calculate overlap % between 2 bed files, need bedtools instal
