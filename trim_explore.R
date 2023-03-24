setwd("/media/barbitoff/DATA/Working issues/BI_Science/variant_calling")
library(ggplot2)
library(ggsci)
library(reshape2)
library(cowplot)
library(RColorBrewer)

trim_data <- read.table("./trimming_r5.tsv", sep='\t', header=T)
trim_data = trim_data[trim_data$Subset == "GRCh37_proteincoding_only.bed.gz", ]
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

# Comparing WES vs. WGS trimming effects
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

plot_grid(trimmers, exptypes, nrow=1, rel_widths = c(1, 0.7), labels=c("a", "b"))

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

### Figure 1 c-d - plotting metric changes against percentage of adapter bases
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
  geom_hline(yintercept = 0, lwd=1, lty=2, col='red') +
  facet_wrap(~variable, nrow=2) + guides(fill=F)
print(mdiff_adapt)


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

caller_recall = aggregate(METRIC.Recall~CallerFilter+Type, caller_data, median)
caller_recall$METRIC.Recall_LE = aggregate(METRIC.Recall~CallerFilter+Type, 
                        caller_le_data, median)$METRIC.Recall
caller_recall$cov_diff = caller_recall$METRIC.Recall - caller_recall$METRIC.Recall_LE
caller_recall

caller_le_data$METRIC.F1_Score_HC = caller_data$METRIC.F1_Score
caller_le_data$F1_diff = caller_le_data$METRIC.F1_Score - 
                         caller_le_data$METRIC.F1_Score_HC
f1_cdiff <- ggplot(caller_le_data, aes(x=CallerFilter, y=F1_diff, fill=CallerFilter)) +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette='Set3') +
  geom_hline(yintercept = 0, lwd=1, lty=2, col='red') + facet_wrap(~Type, ncol=2) + 
  theme(axis.text.x=element_blank(), legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
f1_cdiff

le_f1 <- ggplot(caller_le_data, aes(x=CallerFilter, y=METRIC.F1_Score, 
                                    fill=CallerFilter)) +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette='Set3') +
  facet_wrap(~Type, ncol=2) +  
  theme(axis.text.x=element_blank(), legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
le_f1


####### Paper-ready figures ###################################################

# Figure 1
plot_grid(plot_grid(tbases, treads, nrow=2, labels=c('a', 'b')),
          exptypes, mdiff_adapt, labels=c('', 'c', 'd'), ncol=3,
          rel_widths = c(0.5, 1, 0.7))

# Figure 2
plot_grid(callers, callers_le, 
          plot_grid(f1_cdiff, le_f1, ncol=2, labels=c('c', 'd')),
          nrow=3, labels=c('a', 'b', ''), rel_heights = c(1, 1, 1.2))
