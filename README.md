# Inferring the genetic ancestry of individuals using a labeled training set and principal components analysis

This repostotry contains method to analyze/infer genetic ancestry of indiviuals. First, it utilizes `VCFtools`  to filter the  `vcf file` to high quality variant. Second, it uses `plink2` for linkage disequilibrium prunning and principal components analysis. Last, it applies predictive modeling `random forest` to infer the genetic ancestry of individuals with missing ancestry labels. 

The initial filtering of high quality vairant is compiled inside [`scripts/run_vcf_plink_filter.py`](scripts/run_vcf_plink_filter.py)
The analysis is compiled inside [`scripts/infer_ancestry_label.py`](scripts/infer_ancestry_label.py). 
The pca plots are located inside [`output/figures`](output/figures)
PCA values for the call set [`output/acs_mini_project.filtered.r2.eigenvec`](output/acs_mini_project.filtered.r2.eigenvec)

NOTE:
Most of the output files are not push to this repo due to size. 

Run python scripts on terminal:
Inside the ancestryLabel folder [`python3 scripts/run_vcf_plink_filter.py`] and [`python3 scripts/infer_ancestry_label.py`]

## Tools need to run the analysis
- [VCFtools](https://vcftools.github.io/man_latest.html)
- [bcftools](https://samtools.github.io/bcftools/bcftools.html#head)
- [plink2](https://www.cog-genomics.org/plink/2.0/)
- [vcflib](https://github.com/vcflib/vcflib)
