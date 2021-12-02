
BBTOOLS=/YOURPATH/bbmap_37.23/
BWA=/YOURPATH/bwa
SAMBAMBA=/YOURPATH/sambamba
SAMTOOLS=/YOURPATH/samtools
PICARD=/YOURPATH/picard.jar 
BEDTOOLS=/YOURPATH/bedtools2-2.24.0/bin

##### Preprocessing the reads (each sample separately):
```{bash}

for sample in `cat e02254-20_ISHIV1accessions.txt`   # file with samples' accession numbers
do
  
  NAME=${sample:0:11}
  IN1=${NAME}_1.fastq.gz   # name of file containing downloaded sequence data for a sample
  IN2=${NAME}_2.fastq.gz

  # fix read naming:
  $BBTOOLS/reformat.sh -Xmx60g in=$IN1 in2=$IN2 out=${NAME}_1_fixname.fastq.gz out2=${NAME}_2_fixname.fastq.gz trimreaddescription=t addslash=t slashspace=f
  
  # for reads 1, remove everything left of viral sequences and everything right of linker sequence, then remove ends with quality < 10:
  $BBTOOLS/bbduk2.sh -Xmx60g in=${NAME}_1_fixname.fastq.gz out=${NAME}_1_decontam.fastq.gz lref=NC_001802.1.fasta rref=linker.fasta qtrim=rl trimq=10 threads=10
  # for reads 2, remove everything left of linker and everything right of viral sequence, then remove ends with quality < 10:
  $BBTOOLS/bbduk2.sh -Xmx60g in=${NAME}_2_fixname.fastq.gz out=${NAME}_2_decontam.fastq.gz lref=linker.fasta rref=NC_001802.1.fasta mink=9 qtrim=rl trimq=10 threads=10
  
  # re-pair reads:
  $BBTOOLS/repair.sh -Xmx60g in=${NAME}_1_decontam.fastq.gz in2=${NAME}_2_decontam.fastq.gz out=${NAME}_1_clean.fastq.gz out2=${NAME}_2_clean.fastq.gz repair=t
  
done

```

##### Mapping and filtering (each sample separately):
```{bash}

for sample in `cat e02254-20_ISHIV1accessions.txt`
do
  NAME=${sample:0:11}
  
  # unzipping and changing the names because bwa requires unzipped files with .fq extension:
  gunzip -k ${NAME}_1_clean.fastq.gz
  mv ${NAME}_1_clean.fastq ${NAME}_R1.fq
  gunzip -k ${NAME}_2_clean.fastq.gz
  mv ${NAME}_2_clean.fastq ${NAME}_R2.fq
  
  # map to human genome and filter to retain only reliable mappings:
  $BWA mem -t 10 -M hg19.p13.plusMT.no_alt_analysis_set.fa ${NAME}_R1.fq ${NAME}_R2.fq | $SAMTOOLS view -u -@ 10 - | $SAMTOOLS sort --output-fmt BAM -m 5G -@ 10 - > ${NAME}_hg19.sorted.bam
  $SAMBAMBA view -F "mapping_quality >= 30 and proper_pair and template_length < 900 and not (unmapped or secondary_alignment or supplementary)" -f bam -t 10 ${NAME}_hg19.sorted.bam -o ${NAME}_hg19.filtered.bam
  
  # map to HIV-1 genome and extract reads which pass filters:
  $BWA mem -t 10 -M HIV1 ${NAME}_R1.fq ${NAME}_R2.fq | $SAMTOOLS view -u -@ 10 - | $SAMTOOLS sort --output-fmt BAM -m 5G -@ 10 - > ${NAME}_HIV1.sorted.bam
  $SAMBAMBA view -F "mapping_quality >= 30 and not (unmapped or secondary_alignment or supplementary)" -f bam -t 10 ${NAME}_HIV1.sorted.bam -o ${NAME}_HIV1.filtered.bam
  $SAMTOOLS view ${NAME}_HIV1.filtered.bam | awk '{if(length($10)>20)print $1}' | sort | uniq > ${NAME}_HIV1_reads.txt
  
  # remove reads obtained in previous step from hg19 mapping:
  NRREADS=`wc -l < ${NAME}_HIV1_reads.txt`
  if [[ $NRREADS -gt 0 ]]
    then
      java -jar $PICARD FilterSamReads -INPUT ${NAME}_hg19.filtered.bam -OUTPUT ${NAME}_hg19.nohivreads.bam -READ_LIST_FILE ${NAME}_HIV1_reads.txt -FILTER excludeReadList -TMP_DIR ./tmp
    else
      mv ${NAME}_hg19.filtered.bam ${NAME}_hg19.nohivreads.bam
  fi
  
  # deduplicate:
  $SAMBAMBA markdup -r -t 10 --tmpdir=tmp/ ${NAME}_hg19.nohivreads.bam ${NAME}_hg19.clean.bam
  
done

# keep only read R1 out of all pairs and convert to .bed for downstream analysis:
for bamfile in *_hg19.clean.bam
do
  
  NAME=${bamfile%.bam}
  # get just R1:
  $SAMTOOLS view -@ 1 -f 0x40 -O "BAM" -o ${NAME}_R1.bam ${NAME}.bam
  # turn this to bed:
  $BEDTOOLS/bamToBed -ed -i ${NAME}_R1.bam > ${NAME}_R1.bed

done

```


