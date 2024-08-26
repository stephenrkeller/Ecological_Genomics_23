#!/bin/bash
######
#
# quantify each sample with salmon
#
#######

# -i points to the index files already created
# -l A tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.)
# -p 8 says uses 8 threads
# -o indicates the directory and name of output
# seqbias corrects for random hexamer priming
# gcbias corrects for gcbias, but only when present.

conda activate salmon

for i in $(ls /data/project_data/RNAseq/cleandata | grep '.fq.gz' | cut -f 1-3 -d "_"| uniq);
do

    echo "starting sample ${i}"
    #starting with only name of rep. need to pull out files

    read1=$(ls /data/project_data/RNAseq/cleandata | grep ${i} | grep '_1.qc.fq.gz')
    read2=$(ls /data/project_data/RNAseq/cleandata | grep ${i} | grep '_2.qc.fq.gz')

    salmon quant -i /data/project_data/RNAseq/assembly/hudsonica_index \
        -l A \
         -1 /data/project_data/RNAseq/cleandata/${read1} \
         -2 /data/project_data/RNAseq/cleandata/${read2} \
         -p 6  \
         --softclip \
         --seqBias \
         --gcBias \
         -o /data/project_data/RNAseq/salmon/transcripts_quant_2/${i}

    echo "sample ${i} done"

done
