
#rev project

## Directory structure
The scripts assume the following directories exist:
* work (for intermediate analysis output)
* data (fastq files go here)
* out (for final output)

## Scripts
### starAlign.bash
A script to run STAR aligner on .fastq files.
### checkGenes.R
A simple script to check several housekeeping genes for reasonable output (i.e. no obvious annotation mixups or RNA-Seq breakdowns)
### findIntronExons.R
Analyze gr38 genome to break out regions exclusively exon, intron or mixed and count reads in each region for each sample.
### analyzeIntronExons.R
Calculate differences between treatment and controls in the annotated regions and look for signficant differences.
### functions.R
Helper functions used in the analysis.
