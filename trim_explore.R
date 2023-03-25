setwd("/media/barbitoff/DATA/Working issues/BI_Science/variant_calling")
library(ggplot2)
library(ggsci)
library(reshape2)
library(cowplot)
library(RColorBrewer)

full_data <- read.table("./trimming_r5.tsv", sep='\t', header=T)

### Getting CDS data and pre-processing metric values
trim_data = full_data[full_data$Subset == "GRCh37_proteincoding_only.bed.gz", ]
trim_data = trim_data[trim_data$Filter == "PASS" & trim_data$Subtype == "*", ]
trim_data$METRIC.Precision = as.numeric(trim_data$METRIC.Precision)
trim_data$METRIC.Recall = as.numeric(trim_data$METRIC.Recall)
trim_data$METRIC.F1_Score = as.numeric(trim_data$METRIC.F1_Score)

notrim_data = trim_data[trim_data$Trimming == 'NOTRIM', ]
trim_data = trim_data[trim_data$Trimming != "NOTRIM", ]

rownames(notrim_data) <- apply(notrim_data, 1, 
            function(x) paste(x[1], x[2], x[4], x[5], x[6], sep='_'))

get_diff <- function(x, metric_name, notrim_data=notrim_data) {
  nt_value = as.numeric(notrim_data[paste(x[1], x[2], x[4], x[5], x[6], sep='_'), 
                         metric_name])
  diff = as.numeric(x[metric_name]) - nt_value
  return(diff)
}

trim_data$precision_diff = apply(trim_data, 1, function(x) 
                              get_diff(x, 'METRIC.Precision', notrim_data))
trim_data$recall_diff = apply(trim_data, 1, function(x) 
                              get_diff(x, 'METRIC.Recall', notrim_data))
trim_data$f1_diff = apply(trim_data, 1, function(x) 
                              get_diff(x, 'METRIC.F1_Score', notrim_data))

### Trimmed base count information
base_counts = read.table('base_counts.tsv', sep=' ', header=F)
head(base_counts)
colnames(base_counts) = c('sample', 'trimmer', 'type', 'n_bases')
base_counts = base_counts[base_counts$n_bases > 0, ]

fp_bases = base_counts[base_counts$trimmer == 'FASTP', ]
fp_bases$n_raw = base_counts[base_counts$trimmer == 'NOTRIM', 'n_bases']
fp_bases$trimmed_pct = 1 - fp_bases$n_bases/fp_bases$n_raw
tbases <- ggplot(fp_bases, aes(x=type, y=trimmed_pct, fill=type)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter() + 
  theme_bw() + scale_y_continuous(limits=c(0, 0.07)) +
  guides(fill=F)
print(tbases)

aread_counts = read.table('adapter_reads_counts.txt', sep='\t', header=F)
head(aread_counts)
colnames(aread_counts) = c('sample', 'type', 'trimmed', 'total', 'trimmed_pct')
treads <- ggplot(aread_counts, aes(x=type, y=trimmed_pct, fill=type)) + 
  geom_boxplot() + geom_jitter() + theme_bw() + guides(fill=F)
plot_grid(tbases, treads, nrow=2)


### Comparing WES vs. WGS trimming effects
fastp_data = trim_data[trim_data$Trimming == "FASTP" &
                       trim_data$CallerFilter == 'DV_STANDART' &
                       trim_data$ExperimentType != "LOWEXOME", ]
fd <- melt(fastp_data, measure.vars=c("recall_diff", "precision_diff", "f1_diff"),
           id.vars=c("Sample", "ExperimentType", "Aligner", "CallerFilter", "Type"))
head(fd)
exptypes <- ggplot(fd, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=1) + geom_hline(yintercept = 0, lty=2, col='red') +
  theme_bw() + facet_grid(vars(Type), vars(ExperimentType)) + scale_fill_simpsons() +
  theme(legend.position = 'top', axis.text.x=element_text(angle=45, hjust=1)) +
  guides(fill=F)
print(exptypes)


### Comparing the effects of different trimmer tools
wes_trim = trim_data[trim_data$ExperimentType == 'EXOME' & 
                       trim_data$CallerFilter == 'DV_STANDART', ]
td <- melt(wes_trim, measure.vars=c("recall_diff", "precision_diff", "f1_diff"),
           id.vars=c("Sample", "Trimming", "Aligner", "CallerFilter", "Type"))
head(td)
trimmers <- ggplot(td, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=1) + geom_hline(yintercept = 0, lty=2, col='red') +
  theme_bw() + facet_grid(vars(Type), vars(Trimming)) + scale_fill_simpsons() +
  theme(legend.position = 'top', axis.text.x=element_text(angle=45, hjust=1)) +
  guides(fill=F)
print(trimmers)
plot_grid(trimmers, exptypes, nrow=1, rel_widths = c(1, 0.7), labels=c("a", "b"))


### Plotting metric changes against percentage of adapter bases
fastp_indels = fastp_data[fastp_data$Type == 'INDEL', ]
rownames(fp_bases) = paste(fp_bases$sample, fp_bases$type, sep='_')
rownames(aread_counts) = paste(aread_counts$sample, aread_counts$type, sep='_')
fastp_indels$adapter_bases = sapply(paste(fastp_indels$Sample, 
                                          fastp_indels$ExperimentType,
                                          sep='_'), function(x)
                                    as.numeric(fp_bases[x, 'trimmed_pct']))
fpi = melt(fastp_indels, id.vars=c('Sample', 'ExperimentType', 'adapter_bases'), 
           measure.vars=c('precision_diff', 'recall_diff'))
mdiff_adapt <- ggplot(fpi, aes(x=adapter_bases, y=value, fill=ExperimentType)) + 
  geom_point(size=3, col='black', pch=21) + theme_bw() +
  theme(panel.grid=element_blank()) +
  geom_hline(yintercept = 0, lwd=1, lty=2, col='red') +
  facet_wrap(~variable, nrow=1) + guides(fill=F)
print(mdiff_adapt)


### Testing the effects of trimming in boundaries
pad_data = full_data[grepl('vicinity', full_data$Subset), ]
pad_data = pad_data[pad_data$Subtype == '*' &
                      pad_data$Filter == "PASS" &
                      pad_data$ExperimentType == "EXOME" &
                      pad_data$Caller == "DV_STANDART" &
                      pad_data$Trimming %in% c("FASTP", "NOTRIM"), ]
pad_trim = pad_data[pad_data$Trimming == 'FASTP', ]
pad_notrim = pad_data[pad_data$Trimming == "NOTRIM", ]
rownames(pad_notrim) <- apply(pad_notrim, 1, 
            function(x) paste(x[1], x[6], x[8], sep='_'))
get_diff_pad <- function(x, metric_name, notrim_data=notrim_data) {
  nt_value = as.numeric(notrim_data[paste(x[1], x[6], x[8], sep='_'), 
                                    metric_name])
  diff = as.numeric(x[metric_name]) - nt_value
  return(diff)
}
pad_trim$precision_diff = apply(pad_trim, 1, function(x) 
          get_diff_pad(x, 'METRIC.Precision', pad_notrim))
pad_trim$recall_diff = apply(pad_trim, 1, function(x) 
          get_diff_pad(x, 'METRIC.Recall', pad_notrim))
pad_trim$f1_diff = apply(pad_trim, 1, function(x) 
          get_diff_pad(x, 'METRIC.F1_Score', pad_notrim))
dist_df = data.frame(sets=unique(pad_trim$Subset),
                     distancia=c(0, 100, 125, 150, 25, 50, 75))
pad_trim$Distance = sapply(as.character(pad_trim$Subset),
                     function(x) dist_df[dist_df$sets == x, 'distancia']) 

dist_f1 = aggregate(f1_diff~Type+Distance, pad_trim, median)
dist_f1$se = aggregate(f1_diff~Type+Distance, pad_trim, 
                       function(x) sd(x)/sqrt(length(x)))$f1_diff
distances <- ggplot(dist_f1[dist_f1$Distance <= 100, ], 
        aes(x=Distance, y=f1_diff, ymin=f1_diff-se, ymax=f1_diff+se, col=Type)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2) + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_brewer(palette="Set2") +
  scale_x_continuous(limits=c(-1, 101)) +
  ylab('F1 Score difference') 
print(distances)



### Comparing differnt variant calling pipelines - high-cov EXOME
caller_data = trim_data[trim_data$ExperimentType == 'EXOME' & 
                        trim_data$Trimming == 'FASTP', ]
cd = melt(caller_data, measure.vars=c("recall_diff", "precision_diff", "f1_diff"),
          id.vars=c("Sample", "ExperimentType", "Aligner", "CallerFilter", "Type"))
head(cd)
callers <- ggplot(cd, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=1) + geom_hline(yintercept = 0, lty=2, col='red') +
  theme_bw() + facet_grid(vars(Type), vars(CallerFilter)) + scale_fill_simpsons() +
  theme(legend.position = 'top', axis.text.x=element_text(angle=45, hjust=1)) +
  guides(fill=F)
print(callers)


### Comparing the DeepVariant performance for moderate-coverage WES
coverage_data <- trim_data[trim_data$CallerFilter == "DV_STANDART" &
                             trim_data$Trimming == 'FASTP' &
                             trim_data$ExperimentType != "GENOME", ]
covd = melt(coverage_data, measure.vars=c("recall_diff", "precision_diff", "f1_diff"),
            id.vars=c("Sample", "ExperimentType", "Aligner", "CallerFilter", "Type"))
covplt <- ggplot(covd, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=1) + geom_hline(yintercept = 0, lty=2, col='red') +
  theme_bw() + facet_grid(vars(Type), vars(ExperimentType)) + scale_fill_simpsons() +
  theme(legend.position = 'top', axis.text.x=element_text(angle=45, hjust=1)) +
  guides(fill=F)
print(covplt)

plot_grid(plot_grid(tbases, exptypes, labels=c('a', 'b')), 
          plot_grid(trimmers, covplt, labels=c('c', 'd')),
          callers, labels=c('', '', 'e'), nrow=3)


### Callers in the light of lower coverage data
caller_le_data = trim_data[trim_data$ExperimentType == 'LOWEXOME' & 
                          trim_data$Trimming == 'FASTP', ]
cld = melt(caller_le_data, measure.vars=c("recall_diff", "precision_diff", "f1_diff"),
          id.vars=c("Sample", "ExperimentType", "Aligner", "CallerFilter", "Type"))
head(cld)
callers_le <- ggplot(cld, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=1) + geom_hline(yintercept = 0, lty=2, col='red') +
  theme_bw() + facet_grid(vars(Type), vars(CallerFilter)) + scale_fill_simpsons() +
  theme(legend.position = 'top', axis.text.x=element_text(angle=45, hjust=1)) +
  guides(fill=F)
print(callers_le)

plot_grid(callers, callers_le, nrow=2)

### Brief explortion of recall changes between high-cov and moderate-cov WES
caller_recall = aggregate(METRIC.Recall~CallerFilter+Type, caller_data, median)
caller_recall$METRIC.Recall_LE = aggregate(METRIC.Recall~CallerFilter+Type, 
                        caller_le_data, median)$METRIC.Recall
caller_recall$cov_diff = caller_recall$METRIC.Recall - caller_recall$METRIC.Recall_LE
caller_recall


### Exploring the differences in performance (F1) on high-cov and moderate-cov WES
caller_le_data$METRIC.F1_Score_HC = caller_data$METRIC.F1_Score
caller_le_data$F1_diff = caller_le_data$METRIC.F1_Score - 
                         caller_le_data$METRIC.F1_Score_HC
f1_cdiff <- ggplot(caller_le_data, aes(x=CallerFilter, y=F1_diff, fill=CallerFilter)) +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette='Set3') +
  geom_hline(yintercept = 0, lwd=1, lty=2, col='red') + facet_wrap(~Type, ncol=2) + 
  theme(axis.text.x=element_blank(), legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
f1_cdiff


### And a plot comparing different pipelines' performance in moderate-cov WES
le_f1 <- ggplot(caller_le_data, aes(x=CallerFilter, y=METRIC.F1_Score, 
                                    fill=CallerFilter)) +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette='Set3') +
  facet_wrap(~Type, ncol=2) +  
  theme(axis.text.x=element_blank(), legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
le_f1


####### Paper-ready figures ###################################################

# Figure 1
plot_grid(plot_grid(plot_grid(tbases, treads, nrow=2, labels=c('a', 'b')),
          exptypes, ncol=2, rel_widths=c(0.5, 1), labels=c('', 'c')), 
          plot_grid(mdiff_adapt, distances, nrow=1, labels=c('d', 'e'),
          rel_widths = c(1, 0.75)), nrow=2, rel_heights = c(1, 0.8))

# Figure 2
plot_grid(callers, callers_le, 
          plot_grid(f1_cdiff, le_f1, ncol=2, labels=c('c', 'd')),
          nrow=3, labels=c('a', 'b', ''), rel_heights = c(1, 1, 1.2))
