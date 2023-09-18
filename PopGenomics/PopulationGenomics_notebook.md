# Population Genomics lab notebook

## Author: Stephen Keller

### Affiliation: University of Vermont, Dept. Plant of Plant Biology

### E-mail contact: [srkeller\@uvm.edu](mailto:srkeller@uvm.edu){.email}

### Start Date: 09-11-23

### End Date: TBD

### Project Descriptions: This notebook will document my workflow on the bioinformatics of the Population Genomics section of the Fall 2023 Ecological Genomics course.

# Table of Contents:

-   [Entry 1: 2023-09-11](#id-section1)
-   [Entry 2: 2023-09-13](#id-section2)
-   [Entry 3: 2023-09-18](#id-section3)

------    
<div id='id-section1'/>   


### Entry 1: 2023-09-11.   

-  We reviewed the red spruce study system and the exome capture data
-  We discussed the structure of fastq files (DNA sequence, plus the Qscores)
-  Using the program FastQC, we analyzed the quality of the sequencing runs for one file


------    
<div id='id-section2'/>   


### Entry 2: 2023-09-13.  

-  After discussing the FastQC results, we saw good quality sequence data for most of the read length
-  The initial 5 bp or so had more variable base frequencies, and the very end of the reads had slightly lower Q-scores
-  Based on this, we set up an analysis to trim the reads using the `fastp` program
-  We ran the bash script `fastp.sh` for this
-  We looked at the html files produced by `fastp` and compared pre- and post-trimming -- things looked good!
-  We ended the day setting up our read mapping of the trimmed and cleaned reads using `bwa`


------    
<div id='id-section3'/>   


### Entry 3: 2023-XX-XX.