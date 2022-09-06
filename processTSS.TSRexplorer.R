###Load libraries
library("TSRexploreR", lib.loc = "/home/aem11309/Rlibs/")
library("ChIPseeker", lib.loc = "/home/aem11309/Rlibs/")
library("ggplot2")
library("ggrastr", lib.loc = "/home/aem11309/Rlibs/")
library("ggseqlogo")
library("plyranges", lib.loc = "/home/aem11309/Rlibs/")


##threshold=80
##

###Empty tsrexplorer object to store imported data
exp <- tsr_explorer()
getwd()

###import assembly and annotation
assembly <- "/scratch/aem11309/Working/Maize_Genome/Maize_v5/Zm-B73-REFERENCE-NAM-5.0.fa"
annotation <- "/scratch/aem11309/Working/Maize_Genome/Maize_v5/Zm-B73-REFERENCE-NAM-5.0.gtf"

# Load BAM file
bam_file <- "Andy-MaizeShoot1_AdaptorTrimmed_UMIProc_Aligned.sortedByCoord.out.bam"

##create data frame sample sheet, import TSS
samples <- data.frame(sample_name="MaizeShoot_scRNAMerge", file_1=bam_file, file_2=NA)

exp <- tsr_explorer(sample_sheet=samples,
	genome_assembly=assembly,
	genome_annotation=annotation,
)
exp <- tss_import(exp)

###Import BAM file
exp <- import_bams(exp, paired=TRUE, proper_pair=TRUE, remove_duplicate=TRUE, remove_secondary=TRUE)

###Analyze soft-clipped bases
#softclip_histogram(exp) +
#  theme_bw() +
#  scale_fill_viridis_d()

###G correction
exp <- G_correction(exp)

###TSS aggregation
exp <- tss_aggregate(exp)

###Format TSS
exp <- format_counts(exp, data_type="tss")

###Normalize TSS
exp <- normalize_counts(exp, data_type="tss", method="CPM")

###Annotate TSS
exp <- annotate_features(exp, data_type="tss", feature_type="transcript")

###Explore native threshold
plot_threshold_exploration(exp, max_threshold = 50, samples="MaizeShoot_scRNAMerge", point_size=1) + scale_color_viridis_c()

###Apply Threshold
exp <- apply_threshold(exp, threshold=10, n_samples=1)

###TSS Genomic distribution plot
plot_genomic_distribution(exp, data_type="tss", samples="MaizeShoot_scRNAMerge") +
  scale_fill_viridis_d(direction=-1, name="Annotation")

###TSS Feature detection plot
plot_detected_features(exp, data_type="tss", samples="MaizeShoot_scRNAMerge") +
  scale_fill_viridis_d(direction=-1)

###TSS density plot
plot_density(exp, data_type="tss", samples="MaizeShoot_scRNAMerge")

###TSS Heatmap
plot_heatmap(
  exp, data_type="tss", samples="MaizeShoot_scRNAMerge",
  upstream=500, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

###Cluster TSSs
exp <- tss_clustering(exp, threshold=20, max_distance=20, n_samples=1)

###Associate TSSs with TSRs
exp <- associate_with_tsr(exp)

###Annotate TSRs
exp <- annotate_features(exp, data_type="tsr", upstream=1000, downstream=250, feature_type="transcript")

###Mark dominant TSS per TSR
exp <- mark_dominant(exp, data_type="tss")

###Calculate TSR metrics
exp <- tsr_metrics(exp)

###Annotate TSRs
exp <- annotate_features(exp, data_type="tsr")

###Create genomic distribution plot
plot_genomic_distribution(exp, data_type="tsr", samples="MaizeShoot_scRNAMerge") +
  scale_fill_viridis_d(direction=-1, name="Annotation")

###Export Data
tss_export(
  exp,
  samples = "all",
  file_type = "bedgraph",
  out_dir = NA,
  diff_tss = FALSE,
  sep = "\\t",
  genome_assembly = assembly
)

tss_export(
  exp,
  samples = "all",
  file_type = "bigwig",
  out_dir = NA,
  diff_tss = FALSE,
  sep = "\\t",
  genome_assembly = assembly
)

tsr_export(
  exp,
  samples = "all",
  file_type = "bed",
  out_dir = NA,
  diff_tsr = FALSE,
  sep = "\\t"
)

cts_tss <- get_counts(exp, data_type = "tss", samples = "all")
cts_tsr <- get_counts(exp, data_type = "tsr", samples = "all")
write.table(cts_tss, file = "MaizeShoot_scRNAMerge_TSSCountsTable.txt", sep = ",", col.names = NA)
write.table(cts_tsr, file = "MaizeShoot_scRNAMerge_TSRCountsTable.txt", sep = ",", col.names = NA)
