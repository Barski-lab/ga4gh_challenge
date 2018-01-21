# ChIP-Seq CWL pipeline from BioWardrobe 

This workflow is a CWL version of a Python pipeline from BioWardrobe (Kartashov and Barski, 2015). It starts by using BowTie (Langmead et al., 2009) to perform alignment to a reference genome, resulting in an unsorted SAM file. The SAM file is then sorted and indexed with samtools (Li et al., 2009) to obtain a BAM file and a BAI index. Next MACS2 (Zhang et al., 2008) is used to call peaks and to estimate fragment size. In the last few steps, the coverage by estimated fragments is calculated from the BAM file and is reported in bigWig format. The pipeline also reports statistics, such as read quality, peak number and base frequency, and other troubleshooting information using tools such as [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) and [bamtools](https://github.com/pezmaster31/bamtools).

This workflow (v0.0.2) is being used for GA4GH/DREAM challenge (phase 2).

The pipeline scheme is available in the Wiki.

To download workflow with testing input data use
```
git clone --recursive https://github.com/Barski-lab/ga4gh_challenge.git
```
It will clone submodule [ga4gh_challenge_data](https://github.com/michael-kotliar/ga4gh_challenge_data)
into `data`folder.
Uncompress input FASTQ file using
```
cd ga4gh_challenge/data
./prepare_inputs.sh
```
