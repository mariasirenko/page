# page
post alignment genotyping 

Genotyping SNVs
1. Run cellranger count
2. Run getBaseCalls on the cellranger count output directory
2a. plotnumis() to look at number of UMIs per cell
3. Convert basecalls into wildtype and mutant assignments
4. Merge genotype assignments with seurat object

The main idea behind this pipeline is to take advantage of the flexibility that genotyping fastqs after alignment to the genome offers. By first running the cell ranger pipeline,
1. Cell barcodes are matched to the 10X whitelist
2. UMIs are collapsed
3. Reads are aligned to the genome. The STAR aligner is tolerant of mismatches, sequencing errors, splicing, and frameshifts.

We can genotype only the relevant reads at our position of interest (potentially saving time and resources).
We can find where off-target reads are mapping (which genes, etc) to optimize GOT primers in future experiments.
The function getHileup() uses the python library pysam which is a wrapper around samtools and htslib, which are written in C. We use pysam/samtools to open the BAM file created by cell ranger count. In addition to the common BAM file flags, cell ranger output BAM also has custom flags which we use to retreive the cell barcode and UMI for each read.
One key feature of getHileup() is that the output text file contains the UMIs, cell barcodes, and the allele at the pileup position. Once this step is run, numerous downstream analyses or steps (such as additional UMI collapsing - see below) can be run before proceeding to Step 3. Note that cell ranger count is already matching barcodes to the whitelist using 1 hamming distance. We are using the cell barcode that is corrected by cell ranger (flag: CB). Note, it is possible to use the original uncorrected cell barcode instead (flag: CR). Thus, we are only looking at cell barcodes which cell ranger was able to match to the whitelist. More on that here.

Another consideration is UMI collapsing. We are not doing any UMI collapsing ourselves here. The reason for this is that cell ranger is already doing some UMI collapsing by merging UMIs 1-Hamming Distance apart to the UMI with the higher count. More here. However, as stated above, the output of getHileup() contains all the UMIs per cell so one could write a custom UMI collapsing script to run before proceeding to the next steps of wildtype and mutant assignments. We have not done this already because most of our libraries only have 1-2 UMIs per cell. 

For mutant and wildtype cell assignments, a cell is first assigned a genotype of "NA". Then, we look for reads for the cells which are wildtype and overwrite the assignment as "WT". Next we look for the mutant allele in the base calls and if one is found, we overwrite the assignment as "MT". This means that every cell that contains at least one mutant read is called as mutant without consideration for how many wildtype reads there are in that cell. The reason for this is that we believe the false positive rate (caused by sequencing error) is relatively low compared to the false negative rate (allelic dropout in a heterozygous mutant). We also disregard any reads which are not the wildtype or mutant allele for we believe them to be sequencing errors. Due to the high false negative rate, we also suggest that in downstream analyses, such as in seurat, one consider the 'wildtype' label with a grain of salt.
