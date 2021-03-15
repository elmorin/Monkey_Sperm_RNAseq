library(dplyr)
library(plyr)
library(glycanr)
library(viridis)
library(ggrepel)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(tidyverse)
library(Hmisc)
library(VennDiagram)
library(DESeq2)
library(edgeR)
library(limma)
library(janitor)

#exon <- read.table('/Users/elyse.morin/Documents/Monkey Sperm/MacaM_Rhesus_Genome_Annotation_v7.8.2.gtf')
#exon_simp <- exon[exon$V3=="exon",]
#keep <- c("V1","V4","V5","V6","V7","V16")
#exon_simp2 <- exon_simp[,(names(exon_simp) %in% keep)]
#write.table(exon_simp2, file='/Users/elyse.morin/Documents/Monkey Sperm/MacaM_exon.bed', sep='\t', quote=FALSE, col.names=FALSE, row.names = FALSE)

setwd('/Users/elyse.morin/Documents/Dropbox_backup/Monkey Sperm/coverage text files/')

#coverage historgram
dat_30_fwd <- read.table('coverage_results_30_fwd.txt')
dat_31_fwd <- read.table('coverage_results_31_fwd.txt')
dat_32_fwd <- read.table('coverage_results_32_fwd.txt')
dat_33_fwd <- read.table('coverage_results_33_fwd.txt')
dat_34_fwd <- read.table('coverage_results_34_fwd.txt')
dat_35_fwd <- read.table('coverage_results_35_fwd.txt')
dat_36_fwd <- read.table('coverage_results_36_fwd.txt')

dat_30_rev <- read.table('coverage_results_30_rev.txt')
dat_31_rev <- read.table('coverage_results_31_rev.txt')
dat_32_rev <- read.table('coverage_results_32_rev.txt')
dat_33_rev <- read.table('coverage_results_33_rev.txt')
dat_34_rev <- read.table('coverage_results_34_rev.txt')
dat_35_rev <- read.table('coverage_results_35_rev.txt')
dat_36_rev <- read.table('coverage_results_36_rev.txt')

dat_30 <- read.table('coverage_results_30.txt')
dat_31 <- read.table('coverage_results_31.txt')
dat_32 <- read.table('coverage_results_32.txt')
dat_33 <- read.table('coverage_results_33.txt')
dat_34 <- read.table('coverage_results_34.txt')
dat_35 <- read.table('coverage_results_35.txt')
dat_36 <- read.table('coverage_results_36.txt')

# Preprocessing data
key <- distinct(dat_30_fwd[,5:6])

control_fwd <- list(dat_30_fwd,dat_31_fwd,dat_32_fwd,dat_33_fwd)
k <- c(7,8,9,10)
keeps <- c("V6","V7","V8","V9")
control_sums_fwd <- data.frame()
animal <- 30
for (i in control_fwd){
  i[ , k] <- apply(i[ , k], 2,function(x) as.numeric(as.character(x)))
  simp <- i[,(names(i) %in% keeps)]
  sum <- simp %>%
    group_by(V6) %>%
    summarise_each(funs(sum))
  sum$prop <- sum$V8/sum$V9
  sum$subject <- animal
  animal <- animal+1
  control_sums_fwd <- rbind(control_sums_fwd, sum)
}
#add strand back in
control_sums_fwd_stranded <- left_join(control_sums_fwd, key, by = "V6")

control_rev <- list(dat_30_rev,dat_31_rev,dat_32_rev,dat_33_rev)
k <- c(7,8,9,10)
keeps <- c("V6","V7","V8","V9")
control_sums_rev <- data.frame()
animal <- 30
for (i in control_rev){
  i[ , k] <- apply(i[ , k], 2,function(x) as.numeric(as.character(x)))
  simp <- i[,(names(i) %in% keeps)]
  sum <- simp %>%
    group_by(V6) %>%
    summarise_each(funs(sum))
  sum$prop <- sum$V8/sum$V9
  sum$subject <- animal
  animal <- animal+1
  control_sums_rev <- rbind(control_sums_rev, sum)
}
#add strand back in
control_sums_rev_stranded <- left_join(control_sums_rev, key, by = "V6")

control_coverage_sense <- rbind(control_sums_fwd_stranded[control_sums_fwd_stranded$V5=='+',], control_sums_rev_stranded[control_sums_rev_stranded$V5=='-',])
colnames(control_coverage_sense) <- c("gene","num_reads","bases_cov","total_bases","Coverage","subject","strand")
control_coverage_sense$read_type <- "sense"
control_coverage_sense$group <- "control"
control_coverage_antisense <- rbind(control_sums_fwd_stranded[control_sums_fwd_stranded$V5=='-',], control_sums_rev_stranded[control_sums_rev_stranded$V5=='+',])
colnames(control_coverage_antisense) <- c("gene","num_reads","bases_cov","total_bases","Coverage","subject","strand")
control_coverage_antisense$read_type <- "antisense"
control_coverage_antisense$group <- "control"


malt_fwd <- list(dat_34_fwd,dat_35_fwd,dat_36_fwd)
malt_sums_fwd <- data.frame()
animal <- 34
for (i in malt_fwd){
  i[ , k] <- apply(i[ , k], 2,function(x) as.numeric(as.character(x)))
  simp <- i[,(names(i) %in% keeps)]
  sum <- simp %>%
    group_by(V6) %>%
    summarise_each(funs(sum))
  sum$prop <- sum$V8/sum$V9
  sum$subject <- animal
  animal <- animal+1
  malt_sums_fwd <- rbind(malt_sums_fwd, sum)
}
#add strand back in
malt_sums_fwd_stranded <- left_join(malt_sums_fwd, key, by = "V6")

malt_rev <- list(dat_34_rev,dat_35_rev,dat_36_rev)
malt_sums_rev <- data.frame()
animal <- 34
for (i in malt_rev){
  i[ , k] <- apply(i[ , k], 2,function(x) as.numeric(as.character(x)))
  simp <- i[,(names(i) %in% keeps)]
  sum <- simp %>%
    group_by(V6) %>%
    summarise_each(funs(sum))
  sum$prop <- sum$V8/sum$V9
  sum$subject <- animal
  animal <- animal+1
  malt_sums_rev <- rbind(malt_sums_rev, sum)
}
#add strand back in
malt_sums_rev_stranded <- left_join(malt_sums_rev, key, by = "V6")
malt_coverage_sense <- rbind(malt_sums_fwd_stranded[malt_sums_fwd_stranded$V5=='+',], malt_sums_rev_stranded[malt_sums_rev_stranded$V5=='-',])
colnames(malt_coverage_sense) <- c("gene","num_reads","bases_cov","total_bases","Coverage","subject","strand")
malt_coverage_sense$read_type <- "sense"
malt_coverage_sense$group <- "malt"
malt_coverage_antisense <- rbind(malt_sums_fwd_stranded[malt_sums_fwd_stranded$V5=='-',], malt_sums_rev_stranded[malt_sums_rev_stranded$V5=='+',])
colnames(malt_coverage_antisense) <- c("gene","num_reads","bases_cov","total_bases","Coverage","subject","strand")
malt_coverage_antisense$read_type <- "antisense"
malt_coverage_antisense$group <- "malt"

##combine data
coverage_tot <- rbind(control_coverage_sense,control_coverage_antisense,malt_coverage_antisense, malt_coverage_sense)

write.csv(coverage_tot, 'Full_reads_cov_df.csv', row.names=FALSE)



#only genes with average expression >= 100 reads
avg_exp <- aggregate(coverage_tot$num_reads, list(coverage_tot$gene, coverage_tot$subject), sum)
avg_exp_ <- aggregate(avg_exp$x, list(avg_exp$Group.1), mean)
avg_exp_100 <- avg_exp_[avg_exp_$x>=100,]
subset_genes <- avg_exp_100$Group.1

coverage_tot_100 <- coverage_tot[coverage_tot$gene %in% subset_genes,]

coverage_tot_100_sense <- coverage_tot_100[coverage_tot_100$read_type=="sense",]
coverage_tot_100_anti <- coverage_tot_100[coverage_tot_100$read_type=="antisense",]


#normalizing (median normalized expression (log2))
#sense
coverage_tot_100_sense_avg <- aggregate(coverage_tot_100_sense[,2:5], list(coverage_tot_100_sense$gene, coverage_tot_100_sense$group), mean)
colnames(coverage_tot_100_sense_avg) <- c("gene","group","num_reads","bases_cov","total_bases","Coverage")
num_reads_median <- median(coverage_tot_100_sense_avg$num_reads)
coverage_tot_100_sense_avg$norm <- log2(coverage_tot_100_sense_avg$num_reads/num_reads_median)
coverage_tot_100_sense_avg <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg$norm!='-Inf',]

#antisense
coverage_tot_100_anti_avg <- aggregate(coverage_tot_100_anti[,2:5], list(coverage_tot_100_anti$gene, coverage_tot_100_anti$group), mean)
colnames(coverage_tot_100_anti_avg) <- c("gene","group","num_reads","bases_cov","total_bases","Coverage")
num_reads_median <- median(coverage_tot_100_anti_avg$num_reads)
coverage_tot_100_anti_avg$norm <- log2(coverage_tot_100_anti_avg$num_reads/num_reads_median)
coverage_tot_100_anti_avg <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg$norm!='-Inf',]



######
#Coverage x Expression Plots Sense

coverage_tot_100_sense_control <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg$group=="control",]
coverage_tot_100_sense_maltreated <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg$group=="malt",]

sense_control <- coverage_tot_100_sense_control %>%
  group_by(group) %>%
  mutate(rank = order(order(norm, decreasing=TRUE)))
sense_malt <- coverage_tot_100_sense_maltreated %>%
  group_by(group) %>%
  mutate(rank = order(order(norm, decreasing=TRUE)))
ranked_data_sense <- rbind(sense_control, sense_malt)

ranked_data_sense$sperm_genes <- ifelse(ranked_data_sense$gene=="PRM2" | ranked_data_sense$gene=="PRM1" | ranked_data_sense$gene=="TNP2",ranked_data_sense$norm, NA)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))

ggplot(ranked_data_sense, aes(x=rank, y=norm, color=Coverage)) +
  geom_point(aes(color = Coverage)) +
  #scale_colour_gradient(low = "red", high = "blue")+
  sc+
  ggtitle("Sense Reads") + 
  xlab("Rank in data set ordered by expression (per group)") +
  ylab("Median normalized exression (log2)") +
  labs(fill = "Coverage") +
  theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.title =element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  ylim(-11,13) +
  xlim(-500,7000)+
  annotate('text', x = 3500, y = 3, label = 'Maltreated', fontface='bold', size = 5) +
  annotate('text', x = 1000, y = -0.5, label = 'Control', fontface='bold', size = 5) +
  geom_point(aes(x=rank, y=sperm_genes, color=Coverage), colour="black", shape=21, size=2, stroke=1)+
  geom_label_repel(data = subset(ranked_data_sense[ranked_data_sense$group=="malt",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = 0, nudge_x = 500, fontface='bold', size = 4)+
  geom_label_repel(data = subset(ranked_data_sense[ranked_data_sense$group=="control",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = 0, nudge_x = -1000, fontface='bold', size = 4)

#Coverage x Expression Plots Antiense
coverage_tot_100_anti_control <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg$group=="control",]
coverage_tot_100_anti_maltreated <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg$group=="malt",]

anti_control <- coverage_tot_100_anti_control %>%
  group_by(group) %>%
  mutate(rank = order(order(norm, decreasing=TRUE)))
anti_malt <- coverage_tot_100_anti_maltreated %>%
  group_by(group) %>%
  mutate(rank = order(order(norm, decreasing=TRUE)))
ranked_data_anti <- rbind(anti_control, anti_malt)

ranked_data_anti$sperm_genes <- ifelse(ranked_data_anti$gene=="PRM2" | ranked_data_anti$gene=="PRM1" | ranked_data_anti$gene=="TNP2",ranked_data_anti$norm, NA)
ranked_data_anti$sperm_labels <- ifelse(ranked_data_anti$gene=="PRM2" | ranked_data_anti$gene=="PRM1" | ranked_data_anti$gene=="TNP2",ranked_data_anti$gene, NA)

ggplot(ranked_data_anti, aes(x=rank, y=norm, color=Coverage)) +
  geom_point(aes(color = Coverage)) +
  #scale_colour_gradient(low = "green", high = "blue")+
  sc+
  ggtitle("Antisense Reads") + 
  xlab("Rank in data set ordered by expression (per group)") +
  ylab("Median normalized exression (log2)") +
  labs(fill = "Coverage") +
  theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.title =element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  ylim(-11,13) +
  xlim(-500,7000)+
  annotate('text', x = 4000, y = 2, label = 'Maltreated', fontface='bold', size = 5) +
  annotate('text', x = 3000, y = -1.5, label = 'Control', fontface='bold', size = 5)+
  geom_point(aes(x=rank, y=sperm_genes, color=Coverage), colour="black", shape=21, size=2, stroke=2)+
  geom_label_repel(data = subset(ranked_data_anti[ranked_data_anti$group=="malt",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = 0, nudge_x = 1000, fontface='bold', size = 4)+
  geom_label_repel(data = subset(ranked_data_anti[ranked_data_anti$group=="control",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = -0.5, nudge_x = -1000,fontface='bold', size = 4)


####Plotting top 100 genes, just sense reads

top_100_sense <- ranked_data_sense[ranked_data_sense$rank<=100,]

ggplot(top_100_sense, aes(x=rank, y=norm, color=Coverage)) +
  geom_point(aes(color = Coverage)) +
  #scale_colour_gradient(low = "green", high = "blue")+
  sc+
  ggtitle("Sense Reads: top 100") + 
  xlab("Rank in data set ordered by expression (per group)") +
  ylab("Median normalized exression") +
  labs(fill = "Coverage") +
  theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.title =element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  ylim(2,13) +
  xlim(-5,100)+
  annotate('text', x = 90, y = 6.7, label = 'Maltreated', fontface='bold', size = 5) +
  annotate('text', x = 90, y = 4, label = 'Control', fontface='bold', size = 5) +
  geom_point(aes(x=rank, y=sperm_genes, color=Coverage), colour="black", shape=21, size=2, stroke=2)+
  geom_label_repel(data = subset(top_100_sense[top_100_sense$group=="malt",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = 0.5, nudge_x = 10, fontface='bold', size = 4)+
  geom_label_repel(data = subset(top_100_sense[top_100_sense$group=="control",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = -0.75, nudge_x = -8, fontface='bold', size = 4)

####Plotting top 100 genes, just antisense reads

top_100_antisense <- ranked_data_anti[ranked_data_anti$rank<=100,]

ggplot(top_100_antisense, aes(x=rank, y=norm, color=Coverage)) +
  geom_point(aes(color = Coverage)) +
  #scale_colour_gradient(low = "green", high = "blue")+
  sc+
  ggtitle("Antisense Reads: top 100") + 
  xlab("Rank in data set ordered by expression (per group)") +
  ylab("Median normalized exression") +
  labs(fill = "Coverage") +
  theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.title =element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  ylim(2,13) +
  xlim(-5,100)+
  annotate('text', x = 90, y = 7.5, label = 'Maltreated', fontface='bold', size = 5) +
  annotate('text', x = 90, y = 5, label = 'Control', fontface='bold', size = 5) +
  geom_point(aes(x=rank, y=sperm_genes, color=Coverage), colour="black", shape=21, size=2, stroke=2)+
  geom_label_repel(data = subset(top_100_antisense[top_100_antisense$group=="malt",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = 0.5, nudge_x = 10, fontface='bold', size = 4)+
  geom_label_repel(data = subset(top_100_antisense[top_100_antisense$group=="control",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = -0.75, nudge_x = -8, fontface='bold', size = 4)


######################################
#FIGURE 1

#histograms

control_coverage_sense_genes <- control_coverage_sense[,1:2]
control_coverage_sense_expressed <- control_coverage_sense_genes %>%
  group_by(gene) %>%
  summarise_each(funs(mean))
control_sense_genes <- control_coverage_sense_expressed[control_coverage_sense_expressed$num_reads>10,]
keep <- list(control_sense_genes$gene)
control_coverage_sense_nonzeroexp <- control_coverage_sense[control_coverage_sense$gene %in% keep[[1]],]

control_coverage_anti_genes <- control_coverage_antisense[,1:2]
control_coverage_anti_expressed <- control_coverage_anti_genes %>%
  group_by(gene) %>%
  summarise_each(funs(mean))
control_anti_genes <- control_coverage_anti_expressed[control_coverage_anti_expressed$num_reads>10,]
keep <- list(control_anti_genes$gene)
control_coverage_anti_nonzeroexp <- control_coverage_antisense[control_coverage_antisense$gene %in% keep[[1]],]

malt_coverage_sense_genes <- malt_coverage_sense[,1:2]
malt_coverage_sense_expressed <- malt_coverage_sense_genes %>%
  group_by(gene) %>%
  summarise_each(funs(mean))
malt_sense_genes <- malt_coverage_sense_expressed[malt_coverage_sense_expressed$num_reads>10,]
keep <- list(malt_sense_genes$gene)
malt_coverage_sense_nonzeroexp <- malt_coverage_sense[malt_coverage_sense$gene %in% keep[[1]],]

malt_coverage_anti_genes <- malt_coverage_antisense[,1:2]
malt_coverage_anti_expressed <- malt_coverage_anti_genes %>%
  group_by(gene) %>%
  summarise_each(funs(mean))
malt_anti_genes <- malt_coverage_anti_expressed[malt_coverage_anti_expressed$num_reads>10,]
keep <- list(malt_anti_genes$gene)
malt_coverage_anti_nonzeroexp <- malt_coverage_antisense[malt_coverage_antisense$gene %in% keep[[1]],]



par(mfrow=c(2,2))
hist(control_coverage_sense_nonzeroexp$Coverage, 
     breaks=40,
     main = "Distribution of coverage by sense reads across control males", 
     xlab="Proportion",
     ylim=c(0,10000))
hist(control_coverage_anti_nonzeroexp$Coverage, 
     breaks=40,
     main = "Distribution of coverage by antisense reads across control males", 
     xlab="Proportion",
     ylim=c(0,10000))
hist(malt_coverage_sense_nonzeroexp$Coverage, 
     breaks=40,
     main = "Distribution of coverage by sense reads across maltreated males", 
     xlab="Proportion",
     ylim=c(0,10000))
hist(malt_coverage_anti_nonzeroexp$Coverage, 
     breaks=40,
     main = "Distribution of coverage by antisense reads across maltreated males", 
     xlab="Proportion",
     ylim=c(0,10000))

#######################################################
###Scatter correlation plots for Coverage

#normalizing (median normalized expression (log2))
#sense
coverage_tot_100_sense_avg <- aggregate(coverage_tot_100_sense[,2:5], list(coverage_tot_100_sense$gene, coverage_tot_100_sense$group), mean)
colnames(coverage_tot_100_sense_avg) <- c("gene","group","num_reads","bases_cov","total_bases","Coverage")
coverage_tot_100_sense_avg <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg$num_reads>10,]
num_reads_median <- median(coverage_tot_100_sense_avg$num_reads)
coverage_tot_100_sense_avg$norm <- log2(coverage_tot_100_sense_avg$num_reads/num_reads_median)
coverage_tot_100_sense_avg <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg$norm!='-Inf',]

#antisense
coverage_tot_100_anti_avg <- aggregate(coverage_tot_100_anti[,2:5], list(coverage_tot_100_anti$gene, coverage_tot_100_anti$group), mean)
colnames(coverage_tot_100_anti_avg) <- c("gene","group","num_reads","bases_cov","total_bases","Coverage")
coverage_tot_100_anti_avg <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg$num_reads>10,]
num_reads_median <- median(coverage_tot_100_anti_avg$num_reads)
coverage_tot_100_anti_avg$norm <- log2(coverage_tot_100_anti_avg$num_reads/num_reads_median)
coverage_tot_100_anti_avg <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg$norm!='-Inf',]

#Coverage x Expression Plots Sense

coverage_tot_100_sense_control <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg$group=="control",]
coverage_tot_100_sense_maltreated <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg$group=="malt",]

sense_control <- coverage_tot_100_sense_control %>%
  group_by(group) %>%
  mutate(rank = order(order(norm, decreasing=TRUE)))
sense_malt <- coverage_tot_100_sense_maltreated %>%
  group_by(group) %>%
  mutate(rank = order(order(norm, decreasing=TRUE)))
ranked_data_sense <- rbind(sense_control, sense_malt)

#Coverage x Expression Plots Antisense
coverage_tot_100_anti_control <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg$group=="control",]
coverage_tot_100_anti_maltreated <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg$group=="malt",]

anti_control <- coverage_tot_100_anti_control %>%
  group_by(group) %>%
  mutate(rank = order(order(norm, decreasing=TRUE)))
anti_malt <- coverage_tot_100_anti_maltreated %>%
  group_by(group) %>%
  mutate(rank = order(order(norm, decreasing=TRUE)))
ranked_data_anti <- rbind(anti_control, anti_malt)


ggplot(ranked_data_sense, aes(x=Coverage, y=norm, color=group)) +
  geom_point(aes(color = group, stroke = 0.005)) +
  #scale_colour_gradient(low = "red", high = "blue")+
  #sc+
  ggtitle("Sense Reads") + 
  xlab("Coverage") +
  ylab("Median normalized exression (log2)") +
  labs(fill = "Group") +
  geom_smooth(method = "lm", se = FALSE)+
  theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.title =element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12))+
  theme_minimal()+
  scale_color_manual(values=c("black", "red"))+
  ylim(-6,12)

ggplot(ranked_data_anti, aes(x=Coverage, y=norm, color=group)) +
  geom_point(aes(color = group, stroke = 0.005)) +
  #scale_colour_gradient(low = "red", high = "blue")+
  #sc+
  ggtitle("Antisense Reads") + 
  xlab("Coverage") +
  ylab("Median normalized exression (log2)") +
  labs(fill = "Group") +
  geom_smooth(method = "lm", se = FALSE)+
  theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.title =element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12))+
  theme_minimal()+
  scale_color_manual(values=c("black", "red"))+
  ylim(-6,12)

ranked_data_sense_control <- ranked_data_sense[ranked_data_sense$group=='control',]
ranked_data_sense_malt <- ranked_data_sense[ranked_data_sense$group=='malt',]

cor.test(ranked_data_sense_control$Coverage, ranked_data_sense_control$norm)
cor.test(ranked_data_sense_malt$Coverage, ranked_data_sense_malt$norm)

ranked_data_anti_control <- ranked_data_anti[ranked_data_anti$group=='control',]
ranked_data_anti_malt <- ranked_data_anti[ranked_data_anti$group=='malt',]

cor.test(ranked_data_anti_control$Coverage, ranked_data_anti_control$norm)
cor.test(ranked_data_anti_malt$Coverage, ranked_data_anti_malt$norm)


######################################
#FIGURE 2

# Finding functional buckets for top 100 genes

top_100_sense_control <- top_100_sense[top_100_sense$group=='control',]

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/Expressionvrank_Gogenes_sense.pdf",width = 8, height = 5)

ggplot(top_100_sense_control, aes(x=rank, y=norm)) +
  geom_point(color="#27AAE1") +
  ggtitle("Sense Reads: top 100") + 
  xlab("Rank in data set ordered by expression (per group)") +
  ylab("Median normalized exression") +
  theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.title =element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  theme_minimal()+
  ylim(4,12) +
  xlim(-5,100)+
  geom_point(aes(x=rank, y=sperm_genes, color=Coverage), colour="black", shape=21, size=2, stroke=2)+
  geom_label_repel(data = subset(top_100_sense_control[top_100_sense_control$group=="control",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = -0.75, nudge_x = -8, fontface='bold', size = 4)
dev.off()


top_100_anti_control <- top_100_antisense[top_100_antisense$group=='control',]

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/Expressionvrank_Gogenes_anti.pdf",width = 8, height = 5)

ggplot(top_100_anti_control, aes(x=rank, y=norm)) +
  geom_point(color="#2BB673") +
  ggtitle("Antisense Reads: top 100") + 
  xlab("Rank in data set ordered by expression (per group)") +
  ylab("Median normalized exression") +
  theme(plot.title = element_text(hjust = 0.5, size=18),
        axis.title =element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  theme_minimal()+
  ylim(4,12) +
  xlim(-5,100)+
  geom_point(aes(x=rank, y=sperm_genes, color=Coverage), colour="black", shape=21, size=2, stroke=2)+
  geom_label_repel(data = subset(top_100_anti_control[top_100_anti_control$group=="control",], gene == "PRM2" | gene=="PRM1" | gene=="TNP2"), aes(label=gene), nudge_y = -0.75, nudge_x = -8, fontface='bold', size = 4)
dev.off()

#heatmap top 100 buckets sense
heatmap_sense <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/Buckets_top100/GO_overlap_sense_top_100_heatmap.csv', header=TRUE)
heatmap_sense[is.na(heatmap_sense)] <- 0
colnames(heatmap_sense)[1] <- 'gene'

ranked <- cbind(top_100_sense_control[,1], top_100_sense_control[,8])
heatmap_sense_rank <- left_join(ranked, heatmap_sense, by="gene")

heatmap_sense_rank <- heatmap_sense_rank[order(heatmap_sense_rank$rank), ]
heatmap_sense_rank[is.na(heatmap_sense_rank)] <- 0

heatmap_sense_rank <- subset(heatmap_sense_rank, select = -c(rank))
heatmap_sense_rank <- heatmap_sense_rank[c(1,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2)]

m <- as.matrix(heatmap_sense_rank[,-1])

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/Control_sense_top100_go_heatmap.pdf",width = 30, height = 20)
heatmap(t(m),Colv=NA, Rowv=NA, scale='none', col=c('light gray','#27AAE1'))
dev.off()

#heatmap top 100 buckets antisense
heatmap_anti <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/Buckets_top100/GO_overlap_anti_top_100_heatmap.csv', header=TRUE)
heatmap_anti[is.na(heatmap_anti)] <- 0
colnames(heatmap_anti)[1] <- 'gene'

ranked <- cbind(top_100_anti_control[,1], top_100_anti_control[,8])
heatmap_anti_rank <- left_join(ranked, heatmap_anti, by="gene")

heatmap_anti_rank <- heatmap_anti_rank[order(heatmap_anti_rank$rank), ]
heatmap_anti_rank[is.na(heatmap_anti_rank)] <- 0

heatmap_anti_rank <- subset(heatmap_anti_rank, select = -c(rank))
heatmap_anti_rank <- heatmap_anti_rank[c(1,8,7,6,5,4,3,2)]

m <- as.matrix(heatmap_anti_rank[,-1])

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/Control_anti_top100_go_heatmap.pdf",width = 30, height = 20)
heatmap(t(m),Colv=NA, Rowv=NA, scale='none', col=c('light gray','#2BB673'))
dev.off()


#Venn Diagram Sense

# Generate 3 sets of 200 words
sense_control_set <- top_100_sense$gene[top_100_sense$group=='control']
sense_malt_set <- top_100_sense$gene[top_100_sense$group=='malt']


# Chart
venn.diagram(
  x = list(sense_control_set, sense_malt_set),
  category.names = c(" " , " "),
  filename = '/Users/elyse.morin/Dropbox/Monkey Sperm/Sense_top100_venn_diagram.tif',
  output=TRUE,
  imagetype="tiff",
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  height =500, 
  width = 500,
  cat.cex = 0.5,
  cat.default.pos = "text",
  cat.pos = 0,
  cat.dist = 0.1,
  fontfamily = "serif"
  )

control_sense_vd <- as.data.frame(sense_control_set)
maltreated_sense_vd <- as.data.frame(sense_malt_set)

sense_overlap <- as.data.frame(control_sense_vd$sense_control_set[control_sense_vd$sense_control_set %in% maltreated_sense_vd$sense_malt_set])
sense_nonoverlap_control <- as.data.frame(control_sense_vd$sense_control_set[!control_sense_vd$sense_control_set %in% maltreated_sense_vd$sense_malt_set])
sense_nonoverlap_malt <- as.data.frame(maltreated_sense_vd$sense_malt_set[!maltreated_sense_vd$sense_malt_set %in% control_sense_vd$sense_control_set])

#Venn Diagram Antisense

# Generate 3 sets of 200 words
anti_control_set <- top_100_antisense$gene[top_100_antisense$group=='control']
anti_malt_set <- top_100_antisense$gene[top_100_antisense$group=='malt']

# Chart
venn.diagram(
  x = list(anti_control_set, anti_malt_set),
  category.names = c(" " , " "),
  filename = '/Users/elyse.morin/Dropbox/Monkey Sperm/Antisense_top100_venn_diagram.tif',
  output=TRUE,
  imagetype="tiff",
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  height =500, 
  width = 500,
  cat.cex = 0.5,
  cat.default.pos = "text",
  cat.pos = 0,
  cat.dist = 0.1,
  fontfamily = "serif"
)


control_anti_vd <- as.data.frame(anti_control_set)
maltreated_anti_vd <- as.data.frame(anti_malt_set)

anti_overlap <- as.data.frame(control_anti_vd$anti_control_set[control_anti_vd$anti_control_set %in% maltreated_anti_vd$anti_malt_set])
anti_nonoverlap_control <- as.data.frame(control_anti_vd$anti_control_set[!control_anti_vd$anti_control_set %in% maltreated_anti_vd$anti_malt_set])
anti_nonoverlap_malt <- as.data.frame(maltreated_anti_vd$anti_malt_set[!maltreated_anti_vd$anti_malt_set %in% control_anti_vd$anti_control_set])


#######
#GO analysis

sense_overlap_GO <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/Venn Diagrams/GO_analysis_sense_87_overlapping_genes.csv', header=TRUE)

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/Venn Diagrams/GO_analysis_sense_87_overlapping_genes.pdf",width = 8, height = 8)

ggplot(data=sense_overlap_GO, mapping=aes(x=reorder(GO.biological.process.complete,-log10(P.value)), y=-log10(P.value)))+
  geom_bar(stat='identity')+
  theme_minimal() +
  coord_flip()+
  ggtitle('GO Analysis: 87 Overlapping Sense Genes')+
  theme(axis.title.y = element_blank())
dev.off()

antisense_control_nonoverlap_GO <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/Venn Diagrams/Monkey_antisense_control_nonoverlap_20_GO.csv', header=TRUE)

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/Venn Diagrams/GO_analysis_antisense_20_nonoverlapping_control_genes.pdf",width = 30, height = 8)

ggplot(data=antisense_control_nonoverlap_GO[1:20,], mapping=aes(x=reorder(GO.biological.process.complete,-log10(P.value)), y=-log10(P.value)))+
  geom_bar(stat='identity')+
  theme_minimal() +
  coord_flip()+
  ggtitle('GO Analysis: 20 Non-overlapping Antisense Control Genes (top 20)')+
  theme(axis.title.y = element_blank())
dev.off()












######################################
#FIGURE 3

#DESeq/LimmaVoom

##Select 100% Coverage ; only genes that overlap control and maltreated

#coverage_100_anti <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg $Coverage==1& coverage_tot_100_anti_avg$num_reads>10,]
#coverage_100_anti_genes <- as.data.frame(coverage_100_anti[1])
#coverage_100_anti_genes_dupes <- coverage_100_anti_genes %>% get_dupes(gene)
#anti_dupes <- as.array(coverage_100_anti_genes_dupes$gene)

#coverage_100_sense <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg $Coverage==1& coverage_tot_100_sense_avg$num_reads>10,]
#coverage_100_sense_genes <- as.data.frame(coverage_100_sense[1])
#coverage_100_sense_genes_dupes <- coverage_100_sense_genes %>% get_dupes(gene)
#sense_dupes <- as.array(coverage_100_sense_genes_dupes$gene)

##Select >=90% Coverage ; only genes that overlap control and maltreated

coverage_100_anti <- coverage_tot_100_anti_avg[coverage_tot_100_anti_avg $Coverage>=.8& coverage_tot_100_anti_avg$num_reads>10,]
coverage_100_anti_genes <- as.data.frame(coverage_100_anti[1])
coverage_100_anti_genes_dupes <- coverage_100_anti_genes %>% get_dupes(gene)
anti_dupes <- as.array(coverage_100_anti_genes_dupes$gene)
write.table(anti_dupes, "anti_dupes_80_coverage.txt", row.names = FALSE, quote = FALSE, col.names=FALSE)

coverage_100_sense <- coverage_tot_100_sense_avg[coverage_tot_100_sense_avg $Coverage>=.8& coverage_tot_100_sense_avg$num_reads>10,]
coverage_100_sense_genes <- as.data.frame(coverage_100_sense[1])
coverage_100_sense_genes_dupes <- coverage_100_sense_genes %>% get_dupes(gene)
sense_dupes <- as.array(coverage_100_sense_genes_dupes$gene)
write.table(sense_dupes, "sense_dupes_80_coverage.txt", row.names = FALSE, quote=FALSE, col.names = FALSE)



###limma-voom
setwd('/Users/elyse.morin/Dropbox/Monkey Sperm/LimmaVoom/')

sense_dupes <- as.list(read.table('/Users/elyse.morin/Dropbox/Monkey Sperm/LimmaVoom/sense_dupes_80_coverage.txt'))
anti_dupes <- read.table('/Users/elyse.morin/Dropbox/Monkey Sperm/LimmaVoom/anti_dupes_80_coverage.txt')
anti_dupes <- as.array(anti_dupes$V1)
sense_dupes <- as.array(sense_dupes$V1)


#ct <- read.table("sperm_counts_anti.txt",header=TRUE, check.names = FALSE,row.names=1)
ct <- read.table("sperm_counts_sense.txt",header=TRUE, check.names = FALSE,row.names=1)
ct <- as.matrix(cbind(ct[,6:12]))
cts <- subset(ct, rownames(ct) %in% sense_dupes)

#coldata <- as.matrix(c("control","control","control","control","maltreated","maltreated","maltreated"))
coldata <- as.matrix(c("control","maltreated","control","maltreated","control","maltreated","control"))
colnames(coldata) <- "condition"
rownames(coldata) <- as.matrix(colnames(cts))[,1]


TS <- factor(coldata, levels=c("maltreated", "control"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)

cts <- DGEList(cts)


#filter data
#obtain cpms
cpm <- cpm(cts)
head(cpm)
# Which values in CPM are greater than 0.5?
thresh <- cpm > 0.5
head(thresh)
table(rowSums(thresh))
# We would like to keep genes that have at least 2 TRUES in each row
keep.exprs <- rowSums(thresh)>=2
summary(keep.exprs)

# Subset the rows of cts to keep the more highly expressed genes
cts <- cts[keep.exprs,, keep.lib.sizes=FALSE]
dim(cts)

#look into cpm thresholding
#this normalising takes into account technical bias

# Get log2 counts per million
#Count data is not normally distributed, so if we want to examine the distributions 
#of the raw counts we need to log the counts. Next we’ll use box plots to check the 
#distribution of the read counts on the log2 scale. We can use the cpm function to 
#get log2 counts per million, which are corrected for the different library sizes. 
#The cpm function also adds a small offset to avoid taking log of zero.
logcounts <- cpm(cts,log=TRUE)
# Check distributions of samples using boxplots
#boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
#abline(h=median(logcounts),col="blue")
#title("Boxplots of logCPMs (unnormalised)")


cts <- calcNormFactors(cts, method = "TMM")
str(cts)
dim(cts)
logcounts <- cpm(cts,log=TRUE)
#boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
#abline(h=median(logcounts),col="blue")
#title("Boxplots of logCPMs (normalised)")

#By far, one of the most important plots we make when we analyse 
#RNA-Seq data are MDSplots. An MDSplot is a visualisation of a principle 
#components analysis, which determines the greatest sources of variation 
#in the data. A principle components analysis is an example of an unsupervised 
#analysis, where we don’t need to specify the groups. If your experiment is well 
#controlled and has worked well, what we hope to see is that the greatest sources 
#of variation in the data are the treatments/groups we are interested in.
#plotMDS(cts)

#voom transform
v <- voom(cts,design,plot = TRUE)
#look at the plot

# Fit the linear model
#figure out how to write this using the similar way to a normal linear model

fit <- lmFit(v,design)
names(fit)

#create contrast matrix
cont.matrix <- makeContrasts(CONTROLvsMALT=control - maltreated,levels=design)
cont.matrix
#Now we can apply the contrasts matrix to the fit object to get the statistics and estimated 
#parameters of our comparison that we are interested in. Here we call the contrasts.fit function 
#in limma.
fit.cont <- contrasts.fit(fit, cont.matrix)

#We can use the limma decideTests function to generate a quick summary of DE genes for the contrasts.
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

#The final step is to call the eBayes function, which performs empirical Bayes shrinkage on the 
#variances, and estimates moderated t-statistics and the associated p-values.
fit.cont <- eBayes(fit.cont)

#Plot log residual standard deviation versus average log expression for a fitted microarray linear model.
plotSA(fit.cont, main="Final model: Mean−variance trend")

#The topTable command will always output the top 10 genes by default, even if they are not 
#statistically significant. We can specify the coefficient we are interested in by the name we 
#used in the contrast matrix, or by the column number.
m <- topTable(fit.cont,coef=1,sort.by="P", n="Inf")
mOrdered <- m[order(m$P.Val, decreasing = FALSE),]
head(mOrdered)
write.csv(mOrdered, "Monkey_sperm_DEG_cov800_overlap_sense_null.csv")
mSig <- subset(mOrdered, adj.P.Val < 0.05)
mSig
nrow(mSig)

#To get the full table (i.e. the information for all genes, not just the top 10) we can specify n="Inf"

# We want to highlight the significant genes. We can get this from decideTests.
# And some plots...
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"CONTROLvsMALT"], values = c(-1, 1))
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL)

#MDS plot (similar to PCA?)
lcpm <- cpm(m, log=TRUE)
col.group <- design
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=col.group, col=col.group)
title(main="A. Sample groups")

#Heat plot
library(gplots)
bycond.topgenes <- mOrdered[1:100]
i <- which(v$ %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

##GO analysis 80
sense_overlap_GO_80 <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/LimmaVoom/sense_dupes_80_coverage_GO_output.csv')

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/LimmaVoom/sense_dupes_80_coverage_GO_output.pdf",width = 8, height = 8)

ggplot(data=sense_overlap_GO_80, mapping=aes(x=reorder(GO.biological.process.complete,-log10(P.value)), y=-log10(P.value)))+
  geom_bar(stat='identity')+
  theme_minimal() +
  coord_flip()+
  ggtitle('GO Analysis: Overlapping Sense Genes with >= 80% Coverage')+
  theme(axis.title.y = element_blank())
dev.off()

anti_overlap_GO_80 <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/LimmaVoom/anti_dupes_80_coverage_GO_output.csv')

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/LimmaVoom/anti_dupes_80_coverage_GO_output.pdf",width = 8, height = 8)

ggplot(data=anti_overlap_GO_80, mapping=aes(x=reorder(GO.biological.process.complete,-log10(P.value)), y=-log10(P.value)))+
  geom_bar(stat='identity')+
  theme_minimal() +
  coord_flip()+
  ggtitle('GO Analysis: Overlapping Antisense Genes with >= 80% Coverage')+
  theme(axis.title.y = element_blank())
dev.off()

##GO analysis top 100
sense_overlap_GO_100 <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/Buckets_top100/GO_sense_top_100.csv')

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/Buckets_top100/sense_top100_coverage_GO_output.pdf",width = 8, height = 8)

ggplot(data=sense_overlap_GO_100, mapping=aes(x=reorder(Gene.Set.Name,-log10(p.value)), y=-log10(p.value)))+
  geom_bar(stat='identity')+
  theme_minimal() +
  coord_flip()+
  ggtitle('GO Analysis: Sense Genes top 100')+
  theme(axis.title.y = element_blank())
dev.off()

anti_overlap_GO_100 <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/Buckets_top100/GO_antisense_top_100.csv')

pdf(file = "/Users/elyse.morin/Dropbox/Monkey Sperm/Buckets_top100/anti_top100_GO_output.pdf",width = 8, height = 8)

ggplot(data=anti_overlap_GO_100, mapping=aes(x=reorder(Gene.Set.Name,-log10(p.value)), y=-log10(p.value)))+
  geom_bar(stat='identity')+
  theme_minimal() +
  coord_flip()+
  ggtitle('GO Analysis: Antisense Genes top 100')+
  theme(axis.title.y = element_blank())
dev.off()




######################################
#top 100 sense and antisense, controls only

top_100_antisense_control <- top_100_antisense[top_100_antisense$group=='control',][,1]
top_100_sense_control <- top_100_sense[top_100_sense$group=='control',][,1]

write.csv(top_100_sense_control,'/Users/elyse.morin/Dropbox/Monkey Sperm/top_100_sense_control.csv')
write.csv(top_100_antisense_control,'/Users/elyse.morin/Dropbox/Monkey Sperm/top_100_antisense_control.csv')

top_100_sense_control <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/top_100_sense_control.csv')
top_100_anti_control <- read.csv('/Users/elyse.morin/Dropbox/Monkey Sperm/top_100_antisense_control.csv')



