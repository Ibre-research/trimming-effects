setwd("/media/barbitoff/DATA/Working issues/BI_Science/variant_calling")
library(ggplot2)
library(ggsci)
library(reshape2)
library(cowplot)
library(RColorBrewer)

full_data <- read.table("./trimming_r8.tsv", sep='\t', header=T)

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
  geom_boxplot(outlier.shape = NA) + geom_jitter(pch=21, size=2, col='black') + 
  theme_bw() + scale_y_continuous(limits=c(0, 0.07)) +
  guides(fill=F) + scale_fill_npg() + theme(panel.grid=element_blank())
print(tbases)

aread_counts = read.table('adapter_reads_counts.txt', sep='\t', header=F)
head(aread_counts)
colnames(aread_counts) = c('sample', 'type', 'trimmed', 'total', 'trimmed_pct')
treads <- ggplot(aread_counts, aes(x=type, y=trimmed_pct, fill=type)) + 
  geom_boxplot() + geom_jitter(pch=21, size=2, col='black') + 
  scale_fill_npg() + theme_bw() + guides(fill=F) +
  theme(panel.grid = element_blank())
plot_grid(tbases, treads, nrow=2)


### Comparing WES vs. WGS trimming effects
fastp_data = trim_data[trim_data$Trimming == "FASTP" &
                       trim_data$CallerFilter == 'DV_STANDART' &
                       trim_data$ExperimentType != "LOWEXOME", ]
fd <- melt(fastp_data, measure.vars=c("recall_diff", "precision_diff", "f1_diff"),
           id.vars=c("Sample", "ExperimentType", "Aligner", "CallerFilter", "Type"))
head(fd)
exptypes <- ggplot(fd, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(size=2, pch=21, col='black') + 
  geom_hline(yintercept = 0, lty=2, col='red') +
  theme_bw() + facet_grid(vars(Type), vars(ExperimentType)) + scale_fill_simpsons() +
  theme(legend.position = 'top', axis.ticks.x=element_blank(),
        panel.grid=element_blank(), axis.text.x=element_blank())
print(exptypes)
aggregate(f1_diff~Type+ExperimentType, fastp_data, median)

fastp_alltype = trim_data[trim_data$Trimming == 'FASTP' &
                          trim_data$CallerFilter == 'DV_STANDART', ]
fastp_alltype$tp_diff <- apply(fastp_alltype, 1, function(x) 
  get_diff(x, 'TRUTH.TP', notrim_data))
fastp_agg = aggregate(tp_diff~Sample+ExperimentType, fastp_alltype, sum)
fastp_agg$Type = 'TOTAL'
tp_cm = rbind(fastp_alltype[, c('Sample', 'ExperimentType', 'tp_diff', 'Type')], 
              fastp_agg)
ggplot(tp_cm, aes(x=tp_diff, fill=Type)) + 
  geom_bar(stat='bin', col='black', binwidth=1) + 
  theme_bw() + facet_grid(vars(Type), vars(ExperimentType)) +
  scale_fill_brewer(palette = "Set2") + ylab('Number of samples') +
  xlab('True positive variants (trimmed - untrimmed)') +
  theme(panel.grid=element_blank()) + guides(fill=F) + 
  scale_x_continuous(limits=c(-5, 5))

tp_cnt <- ggplot(tp_cm[tp_cm$ExperimentType != "LOWEXOME", ], 
                 aes(x=tp_diff, fill=Type)) + 
  geom_bar(stat='bin', col='black', binwidth=1) + 
  theme_bw() + facet_grid(vars(ExperimentType), vars(Type)) +
  scale_fill_brewer(palette = "Set2") + ylab('Number of samples') +
  xlab('True positive variants (trimmed - untrimmed)') +
  theme(panel.grid=element_blank()) + guides(fill=F) + 
  scale_x_continuous(limits=c(-5, 5))
tp_cnt

# Explorint the actual numbers
tp_stat = fastp_alltype[, c('Sample', 'ExperimentType', 'Type', 'TRUTH.TP', 'tp_diff',
                         'METRIC.Recall', 'recall_diff')]
tp_stat$notrim_tp = notrim_data[notrim_data$CallerFilter == 'DV_STANDART',
                                'TRUTH.TP']
tp_stat$notrim_rec = notrim_data[notrim_data$CallerFilter == 'DV_STANDART',
                                 'METRIC.Recall']
exome_indel_stat = tp_stat[tp_stat$Type == 'INDEL' &
                           tp_stat$ExperimentType == 'EXOME', ]



### Comparing the effects of different trimmer tools
wes_trim = trim_data[trim_data$ExperimentType == 'EXOME' & 
                       trim_data$CallerFilter == 'DV_STANDART', ]
td <- melt(wes_trim, measure.vars=c("recall_diff", "precision_diff", "f1_diff"),
           id.vars=c("Sample", "Trimming", "Aligner", "CallerFilter", "Type"))
head(td)
trimmers <- ggplot(td, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(size=2, pch=21, col='black') + 
  geom_hline(yintercept = 0, lty=2, col='red') +
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
  geom_point(size=2, col='black', pch=21) + theme_bw() +
  theme(panel.grid=element_blank()) +
  geom_hline(yintercept = 0, lwd=1, lty=2, col='red') +
  facet_wrap(~variable, nrow=1) + guides(fill=F) + scale_fill_npg()
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
distances <- ggplot(dist_f1[dist_f1$Distance <= 100, ], aes(x=Distance, 
        y=f1_diff, ymin=f1_diff-se, ymax=f1_diff+se, fill=Type,  col=Type)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2, pch=21, col='black') + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2") +
  scale_x_continuous(limits=c(-1, 101)) +
  ylab('F1 Score difference') 
print(distances)

dist_with_main <- rbind(dist_f1, c('INDEL', -25, 
      median(fastp_data[fastp_data$ExperimentType == 'EXOME' &
                        fastp_data$Type == 'INDEL', 'f1_diff']), NA),
      c('SNP', -25, median(fastp_data[fastp_data$ExperimentType == 'EXOME' &
                                      fastp_data$Type == 'SNP', 'f1_diff']), NA))
dist_with_main$f1_diff = as.numeric(dist_with_main$f1_diff)
dist_with_main$Distance = as.numeric(dist_with_main$Distance)
dist_with_main$Type = as.factor(dist_with_main$Type)

distances_wcds <- ggplot(dist_with_main[dist_with_main$Distance <= 100, ], 
  aes(x=Distance, y=f1_diff, fill=Type,  col=Type)) +
  geom_line(lwd=1) + #geom_errorbar(width=0.2, lwd=1.1) + 
  geom_point(size=2, pch=21, col='black') + theme_bw() +
  theme(panel.grid=element_blank(), legend.position = "bottom") + 
  scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2") +
  scale_x_continuous(limits=c(-26, 101)) +
  ylab('F1 Score difference') 
print(distances_wcds)



### Comparing differnt variant calling pipelines - high-cov EXOME
caller_data = trim_data[trim_data$ExperimentType == 'EXOME' & 
                        trim_data$Trimming == 'FASTP', ]
cd = melt(caller_data, measure.vars=c("recall_diff", "precision_diff", "f1_diff"),
          id.vars=c("Sample", "ExperimentType", "Aligner", "CallerFilter", "Type"))
head(cd)
callers <- ggplot(cd, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(size=2, pch=21, col='black') + 
  geom_hline(yintercept = 0, lty=2, col='red') +
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
  geom_boxplot(outlier.shape=NA) + geom_jitter(size=2, pch=21, col='black') + 
  geom_hline(yintercept = 0, lty=2, col='red') +
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
  geom_boxplot(outlier.shape=NA, lwd=0.5) + geom_jitter(size=2, pch=21, col='black') + 
  geom_hline(yintercept = 0, lty=2, col='red') +
  theme_bw() + facet_grid(vars(Type), vars(CallerFilter)) + scale_fill_simpsons() +
  theme(legend.position = 'top', axis.text.x=element_text(angle=45, hjust=1)) +
  guides(fill=F)
print(callers_le)

plot_grid(callers, callers_le, nrow=2)

acd = rbind(cd, cld)
head(acd)
callers_full <- ggplot(acd, aes(x=ExperimentType, y=value, fill=CallerFilter)) + 
  geom_boxplot(lwd=0.5) + geom_hline(yintercept = 0, lwd=0.5, lty=2, col='red') +
  theme_bw() + theme(panel.grid=element_blank(), legend.position='top') +
  facet_grid(vars(variable), vars(Type), scales='free') + #coord_flip() +
  scale_y_continuous(limits=c(-0.02, 0.02)) + scale_fill_brewer(palette="Set3")
print(callers_full)

acd_f1 <- acd[acd$variable == 'f1_diff', ]
callers_f1diff <- ggplot(acd_f1, aes(x=ExperimentType, y=value, fill=CallerFilter)) + 
  geom_boxplot(lwd=0.5) + geom_hline(yintercept = 0, lwd=0.5, lty=2, col='red') +
  theme_bw() + theme(panel.grid=element_blank(), legend.position='top') +
  facet_wrap(~Type, ncol=2) + ylab('F1 difference (trimmed - untrimmed)') +
  scale_y_continuous(limits=c(-0.02, 0.04)) + scale_fill_brewer(palette="Set3")
print(callers_f1diff)


#### Separate F1 differences across coverage group
e_f1diff <- ggplot(acd_f1[acd_f1$ExperimentType == 'EXOME', ], 
                   aes(x=CallerFilter, y=value, fill=CallerFilter)) + 
  geom_boxplot(lwd=0.5) + geom_hline(yintercept = 0, lwd=0.5, lty=2, col='red') +
  theme_bw() + theme(panel.grid=element_blank(), legend.position='top', 
        axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~Type, ncol=2) + ylab('F1 difference (trimmed - untrimmed)') +
  scale_y_continuous(limits=c(-0.02, 0.04)) + scale_fill_brewer(palette="Set3") +
  guides(fill=F)
e_f1diff
le_f1diff <- ggplot(acd_f1[acd_f1$ExperimentType == 'LOWEXOME', ], 
                   aes(x=CallerFilter, y=value, fill=CallerFilter)) + 
  geom_boxplot(lwd=0.5) + geom_hline(yintercept = 0, lwd=0.5, lty=2, col='red') +
  theme_bw() + theme(panel.grid=element_blank(), legend.position='top', 
                     axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~Type, ncol=2) + ylab('F1 difference (trimmed - untrimmed)') +
  scale_y_continuous(limits=c(-0.02, 0.04)) + scale_fill_brewer(palette="Set3") +
  guides(fill=F)
le_f1diff


median_changes = aggregate(value~variable+Type+CallerFilter+ExperimentType, 
                           acd, median)
median_changes[median_changes$variable == 'f1_diff', ]

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
aggregate(METRIC.F1_Score~CallerFilter+Type, caller_data, median)
e_f1 <- ggplot(caller_data, aes(x=CallerFilter, y=METRIC.F1_Score, 
                                   fill=CallerFilter)) +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette='Set3') +
  facet_wrap(~Type, ncol=2) +  
  theme(axis.text.x=element_blank(), legend.position = 'bottom',
        panel.grid = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE)) + guides(fill=F)
print(e_f1)

aggregate(METRIC.F1_Score~CallerFilter+Type, caller_le_data, median)
le_f1 <- ggplot(caller_le_data, aes(x=CallerFilter, y=METRIC.F1_Score, 
                                    fill=CallerFilter)) +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette='Set3') +
  facet_wrap(~Type, ncol=2) +  
  theme(axis.text.x=element_blank(), legend.position = 'bottom',
        panel.grid=element_blank(), axis.ticks.x = element_blank()) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE)) + guides(fill=F)
le_f1


all_caller_data <- rbind(caller_data, caller_le_data[, 1:73])
all_f1 <- ggplot(all_caller_data, aes(x=ExperimentType, y=METRIC.F1_Score, 
                                     fill=CallerFilter)) +
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette='Set3') +
  facet_wrap(~Type, ncol=2) + guides(fill=F) + ylab('F1 score (timmed)') +
  theme(panel.grid = element_blank())
all_f1


####### Paper-ready figures ###################################################

# Figure 1
plot_grid(plot_grid(plot_grid(tbases, treads, nrow=2, labels=c('a', '')),
          exptypes, ncol=2, rel_widths=c(0.5, 1), labels=c('', 'b')),
          tp_cnt, 
          plot_grid(mdiff_adapt, distances_wcds, nrow=1, labels=c('d', 'e'),
          rel_widths = c(1, 0.75)), nrow=3, rel_heights = c(1, 0.8, 0.8),
          labels=c('', 'c', ''))

# Figure 2
plot_grid(callers, callers_le, 
          plot_grid(e_f1, le_f1, ncol=2, labels=c('c', 'd')),
          nrow=3, labels=c('a', 'b', ''), rel_heights = c(1, 1, 1.2))

plot_grid(callers_full, 
          plot_grid(e_f1, le_f1, ncol=2, labels=c('b', 'c')),
          nrow=2, labels=c('a', ''), rel_heights = c(1.7, 1))

plot_grid(callers_f1diff, all_f1,
          nrow=2, labels=c('a', 'b'), rel_heights = c(1.4, 1))

plot_grid(e_f1diff, e_f1, le_f1diff, le_f1, labels=c('a', 'b', 'c', 'd'),
          nrow=2)
