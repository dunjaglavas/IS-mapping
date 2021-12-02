
BBTOOLS=/YOURPATH/bbmap_37.23/

##### Filtering based on linker sequence:

# Sample sequence data which includes linker sequence as part of reads can be additionaly filtered based on this. By keeping only reads which have a match to linker sequence, we ensure that randomly sequenced parts are removed and only mappings which most likely represent real integration sites are kept. Next step is performed on samples in which linker sequence was found amongst overrepresented sequences in QC.
```{bash}

for sample in `cat e02254-20_ISHIV1accessions_linkerpresent.txt`
do
  
  NAME=${sample:0:11}
  IN1=${NAME}_1.fastq.gz    # raw reads files (before pre-processing and mapping)
  IN2=${NAME}_2.fastq.gz
  
  # detect only read pairs where R2 has a part of linker sequence present:
  $BBTOOLS/bbduk2.sh -Xmx60g in=$IN1 in2=$IN2 outm=${NAME}1_linkermatch.fastq.gz outm2=${NAME}2_linkermatch.fastq.gz fref=linkers_short_strand.fasta skipr1=t k=17 threads=10
  # get names of those reads:
  zcat ${NAME}1_linkermatch.fastq.gz | grep "@SRR" | cut -d "@" -f 2 | cut -d " " -f 1 | awk '{print $0"/1"}' > ${NAME}_pair_matches_linker.txt
       
done

```

##### Extract discrete areas where reads mapped (called "islands" or "peaks"):
```{r}
library(data.table)
library(stringr)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)

# list all .bed files obtained at the end of procedure in 1_Read_mapping_and_filtering.md and extract sample names:
bedfiles <- list.files(pattern="*.bed$")
samplenames <- substr(bedfiles, 1, 11)

# read in the table which connects sample name with treatment (available in additional_files):
treatment <- fread("treatment.txt", header = FALSE)
names(treatment) <- c("sample", "samplegroup", "conditions")

all_samples_islands <- list()
all_samples_peaks <- data.table()
mean_cov_in_covered_areas <- c()

for(k in 1:length(bedfiles)) {   # process each sample separately
    
    # import .bed file with positions of mapped R1 reads after filtering:
    s <- import(bedfiles[k], format="bed")
    
    # for samples in which filtering by linker sequence match was possible, keep only mappings by R1 reads whose R2 pair contained a match with linker sequence:
    if(samplenames[k] %in% c("SRR12322275", "SRR12322277", "SRR12322278", "SRR12330757", "SRR12330758", "SRR12330759", "SRR12330760", "SRR12330761", "SRR12330762")) {
      linkermatch <- fread(file=paste0(samplenames[k], "_pair_matches_linker.txt"), header=FALSE)
      s <- s[s$name %in% linkermatch$V1]
    }

    # since integration is staggered, starts and ends have to be moved for 2 or 3 bp to match the central bp of integration site:
    s_dt <- as.data.table(s)
    s_dt <- s_dt[, start := ifelse(strand == "+", start+2, start-2)
                 ][, end := ifelse(strand == "+", end+3, end-3)]
    s <- GRanges(s_dt[, score := NULL])
    seqlevels(s, pruning.mode="coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
    seqlengths(s) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
    
    # get coverage, extract areas where there is something mapped (islands), and calculate average coverage on that areas:
    cov <- coverage(s)
    isl <- reduce(s)
    isl$sample <- samplenames[k]
    isl$group <- treatment[sample == samplenames[k]]$samplegroup
    isl$treatment <- treatment[sample == samplenames[k]]$conditions
    mean_cov_in_covered_areas[k] <- sum(sum(cov)) / sum(width(isl))
    
    # find most likely IS position for each island:
    isl_cov <- RleViewsList(rleList = cov, rangesList = isl)
    max_position_left <- unlist(start(viewRangeMaxs(isl_cov)))
    max_position_right <- unlist(end(viewRangeMaxs(isl_cov)))
    isl$IS_position <- ifelse(strand(isl) == "+", max_position_left, max_position_right)
    
    # calculate some coverage stats for each island:
    isl$isl_width <- unlist(width(isl_cov))
    isl$isl_mean_cov <- unlist(viewMeans(isl_cov))
    isl$cov_at_max <- unlist(viewMaxs(isl_cov))
    isl$enrichment <- unlist(isl$isl_mean_cov / mean_cov_in_covered_areas[k])
    
    
    # islands wider than average read length are potentially problematic in downstream analysis (during duplicate removal), IF they are saturated (meaning high coverage throughout) - mark them:
    isl$note <- ifelse(width(isl) > quantile(isl$isl_width)[4], "wide", "narrow")   # mark wide peaks
    highcov_cutoff <- quantile(isl$isl_mean_cov)[4]            # decide on cutoff for "high coverage"
    highcoverage <- GRanges(slice(cov, highcov_cutoff))        # get high coverage areas
    wideandhigh <- subsetByOverlaps(isl[isl$note == "wide"], highcoverage)   # overlap the two
    
    saturated <- lapply(cov[wideandhigh], function(x) {
        high_cov_area <- sum(as.numeric(x) > highcov_cutoff)   # if high-coverage area...
        length_cutoff <- 0.85*length(as.numeric(x))            # ...spans more than 85% of bases in peak...
        high_cov_area >= length_cutoff                         # ...return true
    })
    saturated_ranges <- wideandhigh[unlist(saturated)]         # mark saturated islands
    isl$note[subjectHits(findOverlaps(saturated_ranges, isl))] <- "saturated"
       
    # since IS positions are more or less approximate, we can define them as peaks including +/- 5 around IS:
    peaks <- isl
    ranges(peaks) <- IRanges(start=isl$IS_position, end=isl$IS_position)
    peaks <- trim(resize(peaks, width=11, fix = "center"))
    all_samples_peaks <- rbind(all_samples_peaks, as.data.table(peaks))
    
    # save intermediary data for later use:
    saveRDS(all_samples_peaks, file="HIV1_all_peaks.rds") 
}

```

##### Filter islands/peaks to remove duplicates (those represent sequencing and mapping procedure artifacts, rather than real integration sites):
```{r}
library(data.table)
library(stringr)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)

# read in all peaks: 
all_samples_peaks <- readRDS(file="HIV1_all_peaks.rds")
peaks <- GRanges(all_samples_peaks)
seqlevels(peaks, pruning.mode="coarse") <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlengths(peaks) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)

# bin the genome and mark in which bins peaks/islands are located:
binsize <- 100
set.seed(1234)
bins <- trim(shift(tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), 
                              tilewidth = binsize, cut.last.tile.in.chrom = TRUE), sample(1:(binsize/2), 1)))
binoverlap <- findOverlaps(peaks, bins)
group <- as.data.table(binoverlap)[, .SD[1], by=queryHits]$subjectHits   # takes care of peaks which are in multiple bins

peakindices <- 1:length(peaks)
peakindicesfound <- unique(queryHits(binoverlap))
peaks <- peaks[peakindices %in% peakindicesfound]                        # takes care of peaks which are not in any bins (e.g. at last <100bp of chromosome)

# mark bins which have multiple peaks assigned to them:
peaks$bin <- group
peaks <- as.data.table(peaks)
binindices <- peaks[, unique(bin), by=sample]$V1   # unique(bin) stops it from counting in-sample duplicates
duplicatedbins <- unique(binindices[duplicated(binindices)])

# remove in-sample duplicates:
peaks <- peaks[order(sample, bin, enrichment)][, .SD[.N], by=c("sample", "bin")]

# remove all duplicated peaks if one of them is saturated:
uniquepeaks <- peaks[!(bin %in% duplicatedbins)]   
duplicatedpeaks <- peaks[bin %in% duplicatedbins]
removebins <- duplicatedpeaks[note == "saturated"]$bin
duplicatedpeaks_nosat <- duplicatedpeaks[!(bin %in% removebins)]

# out of all remaining duplicates, keep the one with biggest enrichment:
duplicatedpeaks_clean <- duplicatedpeaks_nosat[order(bin, enrichment)][, .SD[.N], by=bin]

# format the end result into a nice-looking table:
peaks_final <- rbind(uniquepeaks, duplicatedpeaks_clean)[order(sample, bin)]
peaks_final <- 
  peaks_final[, celltype := gsub("_.*", "", gsub("Jurkat_", "", treatment))
              ][, virus := ifelse(celltype == "A77V" | celltype == "N74D", celltype, "WT")
                ][, celltype := ifelse(celltype == "A77V" | celltype == "N74D", "WT", celltype)
                  ][, group := gsub(" ", "_", group)
                    ][, group := gsub("/", "", group)
                      ][, IS_position := start + 5   # center position is considered to be the IS
                        ][, .(sample, treatment, group, celltype, virus, seqnames, IS_position, strand, isl_width, isl_mean_cov, cov_at_max, enrichment)]
fwrite(peaks_final, file="HIV1_IS_clean.txt", sep="\t")   # svae the table

# optionally - if you need the insertion sites in .bed format:
peaks_bed <- fread("HIV1_IS_clean.txt")
peaks_bed <- peaks_bed[, .(seqnames, IS_position, strand, isl_mean_cov, celltype, virus, group)
                       ][strand == "+", ':=' (start = IS_position, end = IS_position+1)
                         ][strand == "-", ':=' (start = IS_position-1, end = IS_position)
                           ][, .(seqnames, start, end, strand, isl_mean_cov, celltype, virus, group)]
fwrite(peaks_bed, "HIV1_IS.bed", sep="\t", col.names = FALSE)

# group treatments from different experiment groups together (e.g. WT from capsid mutant experiments and LEDGFp75 knockout experiments are counted together):
peaks_group <- peaks_bed[, experiment := ifelse(virus == "A77V" | virus == "N74D", "capsid mutant", celltype)
                         ][, experiment := ifelse(celltype == "LKO" | celltype == "LEDGFKO", "LKO", experiment)
                           ][, experiment := ifelse(celltype == "NC" | celltype == "NT", "WT", experiment)
                             ][, experiment := ifelse(celltype == "IBD-/-", "IBD--", experiment)]

# optionally - if you want to split data by experimant conditions in separate .bed files:
peaks_group_list <- split(peaks_group, by="experiment")
for(i in 1:length(peaks_group_list)) {
  peaks_out <- peaks_group_list[[i]][, .(seqnames, start, end, strand, isl_mean_cov)]
  fwrite(peaks_out, str_c("HIV1_", names(h1_list)[i], "_IS.bed"), sep="\t", col.names = FALSE)
}

```
